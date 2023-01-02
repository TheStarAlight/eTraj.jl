using ForwardDiff

"Sample provider which yields electron samples through SFA-AE formula, matching `IonRateMethod=:SFA_AE`"
struct SFAAESampleProvider <: ElectronSampleProvider
    laser           ::Laser;
    target          ::SAEAtom;          # SFA-AE only supports [SAEAtom].
    tSamples        ::AbstractVector;
    ss_pdSamples    ::AbstractVector;
    ss_pzSamples    ::AbstractVector;
    phaseMethod     ::Symbol;           # currently supports :CTMC, :QTMC, :SCTS.
    ionRatePrefix   ::Symbol;           # currently supports :ExpRate.
    function SFAAESampleProvider(;  laser               ::Laser,
                                    target              ::SAEAtom,
                                    sample_tSpan        ::Tuple{<:Real,<:Real},
                                    sample_tSampleNum   ::Int,
                                    simu_phaseMethod    ::Symbol,
                                    rate_ionRatePrefix  ::Symbol,
                                    ss_pdMax            ::Real,
                                    ss_pdNum            ::Int,
                                    ss_pzMax            ::Real,
                                    ss_pzNum            ::Int,
                                    kwargs...   # kwargs are surplus params.
                                    )
        F0 = LaserF0(laser)
        Ip = IonPotential(target)
        γ0 = AngFreq(laser) * sqrt(2Ip) / F0
        # check phase method support.
        if ! (simu_phaseMethod in [:CTMC, :QTMC, :SCTS])
            error("[SFAAESampleProvider] Undefined phase method [$simu_phaseMethod].")
            return
        end
        # check IonRate prefix support.
        if ! (rate_ionRatePrefix in [:ExpRate])
            error("[SFAAESampleProvider] Undefined tunneling rate prefix [$rate_ionRatePrefix].")
            return
        end
        # check Keldysh paramater.
        if γ0 ≥ 1.0
            @warn "[SFAAESampleProvider] Keldysh parameter γ=$γ0, adiabatic (tunneling) condition [γ<<1] unsatisfied."
        end
        # check sampling parameters.
        @assert (sample_tSampleNum>0) "[SFAAESampleProvider] Invalid time sample number $sample_tSampleNum."
        @assert (ss_pdNum>0 && ss_pzNum>0) "[SFAAESampleProvider] Invalid pd/pz sample number $ss_pdNum/$ss_pzNum."
        # finish initialization.
        return new(laser, target,
                range(sample_tSpan[1],sample_tSpan[2];length=sample_tSampleNum),
                range(-abs(ss_pdMax),abs(ss_pdMax);length=ss_pdNum), range(-abs(ss_pzMax),abs(ss_pzMax);length=ss_pzNum),
                simu_phaseMethod, rate_ionRatePrefix
                )
    end
end

"Gets the total number of batches."
function batchNum(sp::SFAAESampleProvider)
    return length(sp.tSamples)
end

"Generates a batch of electrons of `batchId` from `sp` using SFA-AE method."
function generateElectronBatch(sp::SFAAESampleProvider, batchId::Int)
    t = sp.tSamples[batchId]
    Fx::Function = LaserFx(sp.laser)
    Fy::Function = LaserFy(sp.laser)
    dFxt = ForwardDiff.derivative(Fx,t)
    dFyt = ForwardDiff.derivative(Fy,t)
    Fxt = Fx(t)
    Fyt = Fy(t)
    Ft  = hypot(Fxt,Fyt)
    φ   = atan(-Fyt,-Fxt)
    Ip  = IonPotential(sp.target)
    pdNum, pzNum = length(sp.ss_pdSamples), length(sp.ss_pzSamples)
    dim = (sp.phaseMethod == :CTMC) ? 8 : 9 # x,y,z,px,py,pz,t0,rate[,phase]
    rate::Function =
        if  sp.ionRatePrefix == :ExpRate
            ADKRateExp(sp.target)
        else
            #TODO: Add support for other prefixes.
        end
    F2eff(px,py) = Ft^2 - (px*dFxt+py*dFyt)  # F2eff=F²-p⟂⋅F'
    r_exit(px,py) = Ft/2*(px^2+py^2+2Ip)/F2eff(px,py)
    sample_count_thread = zeros(Int,nthreads())
    init_thread = zeros(Float64, dim, nthreads(), pdNum*pzNum) # initial condition (support for multi-threading)

    @threads for ipd in 1:pdNum
        for ipz in 1:pzNum
            pd0, pz0 = sp.ss_pdSamples[ipd], sp.ss_pzSamples[ipz]
            px0 = pd0*-sin(φ)
            py0 = pd0* cos(φ)
            F2eff_ = F2eff(px0,py0)
            if F2eff_ > 0
                sample_count_thread[threadid()] += 1
            else
                continue
            end
            r0 = r_exit(px0,py0)
            x0 = r0*cos(φ)
            y0 = r0*sin(φ)
            z0 = 0.
            init_thread[1:8,threadid(),sample_count_thread[threadid()]] = [x0,y0,z0,px0,py0,pz0,t,rate(sqrt(F2eff_),φ,pd0,pz0)]
            if sp.phaseMethod != :CTMC
                init_thread[9,threadid(),sample_count_thread[threadid()]] = 0.0
            end
        end
    end
    # collect electron samples from different threads.
    init = zeros(Float64, dim, sum(sample_count_thread))
    N = 0
    for i in 1:nthreads()
        init[:,N+1:N+sample_count_thread[i]] = init_thread[:,i,1:sample_count_thread[i]]
        N += sample_count_thread[i]
    end
    return init
end