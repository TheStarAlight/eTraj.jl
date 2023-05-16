
using LinearAlgebra

"Sample provider which yields electron samples through ADK rate formula, matching `IonRateMethod=:ADK`."
struct ADKSampler <: ElectronSampleProvider
    laser           ::Laser;
    target          ::SAEAtomBase;          # ADK only supports [SAEAtomBase].
    monteCarlo      ::Bool;
    tSamples        ::AbstractVector;
    ss_kdSamples    ::AbstractVector;
    ss_kzSamples    ::AbstractVector;
    mc_tBatchSize   ::Integer;
    mc_ktMax        ::Real;
    phaseMethod     ::Symbol;           # currently supports :CTMC, :QTMC, :SCTS.
    ionRatePrefix   ::Symbol;           # currently supports :ExpRate.
    ADKTunExit      ::Symbol;           # currently supports :IpF, :FDM, :Para.
    function ADKSampler(;   laser               ::Laser,
                            target              ::SAEAtomBase,
                            sample_tSpan        ::Tuple{<:Real,<:Real},
                            sample_tSampleNum   ::Int,
                            rate_monteCarlo     ::Bool,
                            simu_phaseMethod    ::Symbol,
                            rate_ionRatePrefix  ::Symbol,
                            adk_ADKTunExit      ::Symbol,
                                #* for step-sampling (!rate_monteCarlo)
                            ss_kdMax            ::Real,
                            ss_kdNum            ::Int,
                            ss_kzMax            ::Real,
                            ss_kzNum            ::Int,
                                #* for Monte-Carlo-sampling (rate_monteCarlo)
                            mc_tBatchSize       ::Int,
                            mc_ktMax            ::Real,
                            kwargs...   # kwargs are surplus params.
                            )
        F0 = LaserF0(laser)
        Ip = IonPotential(target)
        γ0 = AngFreq(laser) * sqrt(2Ip) / F0
        # check phase method support.
        if ! (simu_phaseMethod in [:CTMC, :QTMC, :SCTS])
            error("[ADKSampler] Undefined phase method [$simu_phaseMethod].")
            return
        end
        # check tunneling exit support.
        if ! (adk_ADKTunExit in [:IpF, :FDM, :Para])
            error("[ADKSampler] Undefined tunneling exit method [$adk_ADKTunExit].")
            return
        end
        # switch tunneling exit method if inappropriate.
        if adk_ADKTunExit == :FDM && Ip^2 < 4F0
            adk_ADKTunExit == :IpF
            @warn "[ADKSampler] Tunneling exit method has changed from [FDM] to [IpF] due to failure in tunneling exit determination of [FDM] model."
        elseif adk_ADKTunExit == :Para && Ip^2 < 4*(1-sqrt(Ip/2))*F0
            adk_ADKTunExit == :IpF
            @warn "[ADKSampler] Tunneling exit method has changed from [Para] to [IpF] due to failure in tunneling exit determination of [Para] model."
        end
        # check IonRate prefix support.
        if ! (rate_ionRatePrefix in [:ExpRate, :ExpPre, :ExpJac, :Full])
            error("[ADKSampler] Undefined tunneling rate prefix [$rate_ionRatePrefix].")
            return
        end
        # check Keldysh parameter.
        if γ0 ≥ 0.5
            @warn "[ADKSampler] Keldysh parameter γ=$γ0, adiabatic (tunneling) condition [γ<<1] not sufficiently satisfied."
        elseif γ0 ≥ 1.0
            @warn "[ADKSampler] Keldysh parameter γ=$γ0, adiabatic (tunneling) condition [γ<<1] unsatisfied."
        end
        # check sampling parameters.
        @assert (sample_tSampleNum>0) "[ADKSampler] Invalid time sample number $sample_tSampleNum."
        if ! rate_monteCarlo    # check SS sampling parameters.
            @assert (ss_kdNum>0 && ss_kzNum>0) "[ADKSampler] Invalid kd/kz sample number $ss_kdNum/$ss_kzNum."
        else                    # check MC sampling parameters.
            @assert (sample_tSpan[1] < sample_tSpan[2]) "[ADKSampler] Invalid sampling time span $sample_tSpan."
            @assert (mc_tBatchSize>0) "[ADKSampler] Invalid batch size $mc_tBatchSize."
            @assert (mc_ktMax>0) "[ADKSampler] Invalid sampling ptmax $mc_ktMax."
        end
        # finish initialization.
        return if ! rate_monteCarlo
            new(laser,target,rate_monteCarlo,
                range(sample_tSpan[1],sample_tSpan[2];length=sample_tSampleNum),
                range(-abs(ss_kdMax),abs(ss_kdMax);length=ss_kdNum), range(-abs(ss_kzMax),abs(ss_kzMax);length=ss_kzNum),
                0,0,    # for MC params. pass meaningless values
                simu_phaseMethod,rate_ionRatePrefix,adk_ADKTunExit)
        else
            tSamples = rand(sample_tSampleNum) .* (sample_tSpan[2]-sample_tSpan[1]) .+ sample_tSpan[1]
            new(laser,target,rate_monteCarlo,
                tSamples,
                0:0,0:0,    # for SS params. pass meaningless values
                mc_tBatchSize, mc_ktMax,
                simu_phaseMethod,rate_ionRatePrefix,adk_ADKTunExit)
        end
    end
end
"Gets the total number of batches."
function batchNum(sp::ADKSampler)
    return length(sp.tSamples)
end
"Generates a batch of electrons of `batchId` from `sp` using ADK method."
function generateElectronBatch(sp::ADKSampler, batchId::Int)
    t = sp.tSamples[batchId]
    Fx::Function = LaserFx(sp.laser)
    Fy::Function = LaserFy(sp.laser)
    Fxt = Fx(t)
    Fyt = Fy(t)
    Ft = hypot(Fxt,Fyt)
    Ax::Function = LaserAx(sp.laser)
    Ay::Function = LaserAy(sp.laser)
    φ  = atan(-Fyt,-Fxt)
    Z  = AsympNuclCharge(sp.target)
    Ip = IonPotential(sp.target)
    if Ft == 0
        return nothing
    end
    # determining tunneling exit position.
    r_exit =
        if      sp.ADKTunExit == :IpF
            Ip / Ft
        elseif  sp.ADKTunExit == :FDM
            (Ip + sqrt(Ip^2 - 4Ft)) / 2Ft
        elseif  sp.ADKTunExit == :Para
            (Ip + sqrt(Ip^2 - 4*(1-sqrt(Ip/2))*Ft)) / 2Ft
        end
    # determining ADK rate. ADKRate(F,φ,kd,kz)
    ADKRate::Function =
        if      sp.ionRatePrefix == :ExpRate
            ADKRateExp(sp.target)
        else
            ADKRate_Exp = ADKRateExp(sp.target)
            α = 1.0+Z/sqrt(2Ip)
            if  sp.ionRatePrefix == :ExpPre
                function ADKRate_ExpPre(F,φ,kd,kz)
                    return ADKRate_Exp(F,φ,kd,kz) / ((kd^2+kz^2+2Ip)*Ft^2)^(0.5α)
                end
            elseif  sp.ionRatePrefix == :ExpJac
                function ADKRate_ExpJac(F,φ,kd,kz)
                    transform((tr,kd)) = SVector(kd*Fy(tr)/sqrt(Fx(tr)^2+Fy(tr)^2) - Ax(tr), -kd*Fx(tr)/sqrt(Fx(tr)^2+Fy(tr)^2) - Ay(tr))
                    jac = abs(det(ForwardDiff.jacobian(transform,SVector(t,kd))))
                    return ADKRate_Exp(F,φ,kd,kz) * jac
                end
            elseif  sp.ionRatePrefix == :Full
                function ADKRate_Full(F,φ,kd,kz)
                    transform((tr,kd)) = SVector(kd*Fy(tr)/sqrt(Fx(tr)^2+Fy(tr)^2) - Ax(tr), -kd*Fx(tr)/sqrt(Fx(tr)^2+Fy(tr)^2) - Ay(tr))
                    jac = abs(det(ForwardDiff.jacobian(transform,SVector(t,kd))))
                    return ADKRate_Exp(F,φ,kd,kz) * jac / ((kd^2+kz^2+2Ip)*Ft^2)^(0.5α)
                end
            end
        end
    dim = (sp.phaseMethod == :CTMC) ? 8 : 9 # x,y,z,kx,ky,kz,t0,rate[,phase]
    if ! sp.monteCarlo
        kdNum, kzNum = length(sp.ss_kdSamples), length(sp.ss_kzSamples)
        init = zeros(Float64, dim, kdNum, kzNum) # initial condition
        x0 = r_exit*cos(φ)
        y0 = r_exit*sin(φ)
        z0 = 0.
        @threads for ikd in 1:kdNum
            kd0 = sp.ss_kdSamples[ikd]
            kx0 = kd0*-sin(φ)
            ky0 = kd0* cos(φ)
            for ikz in 1:kzNum
                kz0 = sp.ss_kzSamples[ikz]
                init[1:8,ikd,ikz] = [x0,y0,z0,kx0,ky0,kz0,t,ADKRate(Ft,φ,kd0,kz0)]
                if sp.phaseMethod != :CTMC
                    init[9,ikd,ikz] = 0
                end
            end
        end
        return reshape(init,dim,:)
    else
        init = zeros(Float64, dim, sp.mc_tBatchSize)
        ktMax = sp.mc_ktMax
        x0 = r_exit*cos(φ)
        y0 = r_exit*sin(φ)
        z0 = 0.
        @threads for i in 1:sp.mc_tBatchSize
            # generates random (kd0,kz0) inside circle kd0^2+kz0^2=ktMax^2.
            kd0, kz0 = (rand()-0.5)*2ktMax, (rand()-0.5)*2ktMax
            while kd0^2+kz0^2 > ktMax^2
                kd0, kz0 = (rand()-0.5)*2ktMax, (rand()-0.5)*2ktMax
            end
            kx0 = kd0*-sin(φ)
            ky0 = kd0* cos(φ)
            init[1:8,i] = [x0,y0,z0,kx0,ky0,kz0,t,ADKRate(Ft,φ,kd0,kz0)]
            if sp.phaseMethod != :CTMC
                init[9,i] = 0
            end
        end
        return init
    end
end