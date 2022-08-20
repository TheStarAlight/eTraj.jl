
"Sample provider which yields electron samples through ADK rate formula, matching `IonRateMethod=:ADK`."
struct ADKSampleProvider <: ElectronSampleProvider
    laser           ::Laser;
    target          ::SAEAtom;          # ADK only supports [SAEAtom].
    monteCarlo      ::Bool;
    tSamples        ::AbstractVector;
    ss_pdSamples    ::AbstractVector;
    ss_pzSamples    ::AbstractVector;
    mc_tBatchSize   ::Int;
    mc_ptMax        ::Real;
    phaseMethod     ::Symbol;           # currently supports :CTMC, :QTMC, :SCTS.
    ionRatePrefix   ::Symbol;           # currently supports :ExpRate.
    ADKTunExit      ::Symbol;           # currently supports :IpF, :FDM, :ADK.
    function ADKSampleProvider(;laser               ::Laser,
                                target              ::SAEAtom,
                                sample_tSpan        ::Tuple{<:Real,<:Real},
                                sample_tSampleNum   ::Int,
                                rate_monteCarlo     ::Bool,
                                simu_phaseMethod    ::Symbol,
                                rate_ionRatePrefix  ::Symbol,
                                adk_ADKTunExit      ::Symbol,
                                    #* for step-sampling (!rate_monteCarlo)
                                ss_pdMax            ::Real,
                                ss_pdNum            ::Int,
                                ss_pzMax            ::Real,
                                ss_pzNum            ::Int,
                                    #* for Monte-Carlo-sampling (rate_monteCarlo)
                                mc_tBatchSize       ::Int,
                                mc_ptMax            ::Real,
                                kwargs...   # kwargs are surplus params.
                                )
        F0 = LaserF0(laser)
        Ip = IonPotential(target)
        γ0 = AngFreq(laser) * sqrt(2Ip) / F0
        # check phase method support.
        if ! (simu_phaseMethod in [:CTMC, :QTMC, :SCTS])
            error("Undefined phase method [$simu_phaseMethod].")
            return
        end
        # check tunneling exit support.
        if ! (adk_ADKTunExit in [:IpF, :FDM, :Para])
            error("Undefined tunneling exit method [$adk_ADKTunExit].")
            return
        end
        # switch tunneling exit method if inappropriate.
        if adk_ADKTunExit == :FDM && Ip^2 < 4F0
            adk_ADKTunExit == :IpF
            @warn "Tunneling exit method has changed from [FDM] to [IpF] due to failure in tunneling exit determination of [FDM] model."
        elseif adk_ADKTunExit == :Para && Ip^2 < 4*(1-sqrt(Ip/2))*F0
            adk_ADKTunExit == :IpF
            @warn "Tunneling exit method has changed from [Para] to [IpF] due to failure in tunneling exit determination of [Para] model."
        end
        # check IonRate prefix support.
        if ! (rate_ionRatePrefix in [:ExpRate])
            error("Undefined tunneling rate prefix [$rate_ionRatePrefix].")
            return
        end
        # check Keldysh parameter.
        if γ0 ≥ 0.5
            @warn "Keldysh parameter γ=$γ0, adiabatic (tunneling) condition [γ<<1] not sufficiently satisfied."
        elseif γ0 ≥ 1.0
            @warn "Keldysh parameter γ=$γ0, adiabatic (tunneling) condition [γ<<1] unsatisfied."
        end
        # check sampling parameters.
        @assert (sample_tSampleNum>0) "Invalid time sample number $sample_tSampleNum."
        if ! rate_monteCarlo    # check SS sampling parameters.
            @assert (ss_pdNum>0 && ss_pzNum>0) "Invalid pd/pz sample number $ss_pdNum/$ss_pzNum."
        else                    # check MC sampling parameters.
            @assert (sample_tSpan[1] < sample_tSpan[2]) "Invalid sampling time span $sample_tSpan."
            @assert (mc_tBatchSize>0) "Invalid batch size $mc_tBatchSize."
            @assert (mc_ptMax>0) "Invalid sampling ptmax $mc_ptMax."
        end
        # finish initialization.
        return if ! rate_monteCarlo
            new(laser,target,rate_monteCarlo,
                range(sample_tSpan[1],sample_tSpan[2];length=sample_tSampleNum),
                range(-abs(ss_pdMax),abs(ss_pdMax);length=ss_pdNum), range(-abs(ss_pzMax),abs(ss_pzMax);length=ss_pzNum),
                0,0,    # for MC params. pass meaningless values
                simu_phaseMethod,rate_ionRatePrefix,adk_ADKTunExit)
        else
            tSamples = rand(sample_tSampleNum) .* (sample_tSpan[2]-sample_tSpan[1]) .+ sample_tSpan[1]
            new(laser,target,rate_monteCarlo,
                tSamples,
                0:0,0:0,    # for SS params. pass meaningless values
                mc_tBatchSize, mc_ptMax,
                simu_phaseMethod,rate_ionRatePrefix,adk_ADKTunExit)
        end
    end
end
"Gets the total number of batches."
function batchNum(sp::ADKSampleProvider)
    return length(sp.tSamples)
end
"Generates a batch of electrons of `batchId` from `sp` using ADK method."
function generateElectronBatch(sp::ADKSampleProvider, batchId::Int)
    t = sp.tSamples[batchId]
    Fxt::Function = LaserFx(sp.laser)
    Fyt::Function = LaserFy(sp.laser)
    Fx = Fxt(t)
    Fy = Fyt(t)
    F  = hypot(Fx,Fy)
    φ  = atan(-Fy,-Fx)
    Ip = IonPotential(sp.target)
    # determining tunneling exit position.
    r_exit =
        if      sp.ADKTunExit == :IpF
            Ip / F
        elseif  sp.ADKTunExit == :FDM
            (Ip + sqrt(Ip^2 - 4F)) / 2F
        elseif  sp.ADKTunExit == :Para
            (Ip + sqrt(Ip^2 - 4*(1-sqrt(Ip/2))*F)) / 2F
        end
    # determining ADK rate. ADKRate(F,φ,pd,pz)
    ADKRate::Function =
        if      sp.ionRatePrefix == :ExpRate
            ADKRateExp(sp.target)
        else
            #TODO: Add support for other prefixes.
        end
    dim = (sp.phaseMethod == :CTMC) ? 8 : 9 # x,y,z,px,py,pz,t0,rate[,phase]
    if ! sp.monteCarlo
        pdNum, pzNum = length(sp.ss_pdSamples), length(sp.ss_pzSamples)
        init = zeros(Float64, dim, pdNum, pzNum) # initial condition
        x0 = r_exit*cos(φ)
        y0 = r_exit*sin(φ)
        z0 = 0.
        @threads for ipd in 1:pdNum
            pd0 = sp.ss_pdSamples[ipd]
            px0 = pd0*-sin(φ)
            py0 = pd0* cos(φ)
            for ipz in 1:pzNum
                pz0 = sp.ss_pzSamples[ipz]
                init[1:8,ipd,ipz] = [x0,y0,z0,px0,py0,pz0,t,ADKRate(F,φ,pd0,pz0)]
                if sp.phaseMethod != :CTMC
                    init[9,ipd,ipz] = 0
                end
            end
        end
        return reshape(init,dim,:)
    else
        init = zeros(Float64, dim, sp.mc_tBatchSize)
        ptMax = sp.mc_ptMax
        x0 = r_exit*cos(φ)
        y0 = r_exit*sin(φ)
        z0 = 0.
        @threads for i in 1:sp.mc_tBatchSize
            # generates random (pd0,pz0) inside circle pd0^2+pz0^2=ptMax^2.
            pd0, pz0 = (rand()-0.5)*2ptMax, (rand()-0.5)*2ptMax
            while pd0^2+pz0^2 > ptMax^2
                pd0, pz0 = (rand()-0.5)*2ptMax, (rand()-0.5)*2ptMax
            end
            px0 = pd0*-sin(φ)
            py0 = pd0* cos(φ)
            init[1:8,i] = [x0,y0,z0,px0,py0,pz0,t,ADKRate(F,φ,pd0,pz0)]
            if sp.phaseMethod != :CTMC
                init[9,i] = 0
            end
        end
        return init
    end
end