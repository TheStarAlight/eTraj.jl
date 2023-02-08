using StaticArrays
using NLsolve
using QuadGK
using Base.Threads

"Sample provider which yields electron samples through SFA formula, matching `IonRateMethod=:SFA`"
struct SFASampler <: ElectronSampleProvider
    laser           ::Laser;
    target          ::SAEAtomBase;          # SFA only supports [SAEAtomBase].
    tSamples        ::AbstractVector;
    ss_pdSamples    ::AbstractVector;
    ss_pzSamples    ::AbstractVector;
    phaseMethod     ::Symbol;           # currently supports :CTMC, :QTMC, :SCTS.
    ionRatePrefix   ::Symbol;           # currently supports :ExpRate.
    function SFASampler(;   laser               ::Laser,
                            target              ::SAEAtomBase,
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
        # check phase method support.
        if ! (simu_phaseMethod in [:CTMC, :QTMC, :SCTS])
            error("[SFASampleProvider] Undefined phase method [$simu_phaseMethod].")
            return
        end
        # check IonRate prefix support.
        if ! (rate_ionRatePrefix in [:ExpRate])
            error("[SFASampleProvider] Undefined tunneling rate prefix [$rate_ionRatePrefix].")
            return
        end
        # check sampling parameters.
        @assert (sample_tSampleNum>0) "[SFASampleProvider] Invalid time sample number $sample_tSampleNum."
        @assert (ss_pdNum>0 && ss_pzNum>0) "[SFASampleProvider] Invalid pd/pz sample number $ss_pdNum/$ss_pzNum."
        # finish initialization.
        return new(laser, target,
                range(sample_tSpan[1],sample_tSpan[2];length=sample_tSampleNum),
                range(-abs(ss_pdMax),abs(ss_pdMax);length=ss_pdNum), range(-abs(ss_pzMax),abs(ss_pzMax);length=ss_pzNum),
                simu_phaseMethod, rate_ionRatePrefix
                )
    end
end

"Gets the total number of batches."
function batchNum(sp::SFASampler)
    return length(sp.tSamples)
end

"Generates a batch of electrons of `batchId` from `sp` using SFA method."
function generateElectronBatch(sp::SFASampler, batchId::Int)
    tr   = sp.tSamples[batchId]  # real part of time
    Ax::Function = LaserAx(sp.laser)
    Ay::Function = LaserAy(sp.laser)
    Fx::Function = LaserFx(sp.laser)
    Fy::Function = LaserFy(sp.laser)
    Fxtr = Fx(tr)
    Fytr = Fy(tr)
    Ftr  = hypot(Fxtr,Fytr)
    φ   = atan(-Fytr,-Fxtr)
    ω   = AngFreq(sp.laser)
    Ip  = IonPotential(sp.target)
    if Ftr == 0
        return nothing
    end
    pdNum, pzNum = length(sp.ss_pdSamples), length(sp.ss_pzSamples)
    dim = (sp.phaseMethod == :CTMC) ? 8 : 9 # x,y,z,px,py,pz,t0,rate[,phase]
    sample_count_thread = zeros(Int,nthreads())
    init_thread = zeros(Float64, dim, nthreads(), pdNum*pzNum) # initial condition (support for multi-threading)

    "Saddle-point equation."
    function saddlePointEquation(AxF,AyF,tr,ti,pd,pz)
        Axt = AxF(tr+1im*ti)
        Ayt = AyF(tr+1im*ti)
        px =  pd*imag(Ayt)/sqrt(imag(Axt)^2+imag(Ayt)^2) - real(Axt)
        py = -pd*imag(Axt)/sqrt(imag(Axt)^2+imag(Ayt)^2) - real(Ayt)
        return SVector(real(((px+Axt)^2+(py+Ayt)^2+pz^2)/2 + Ip))
    end

    @threads for ipd in 1:pdNum
        for ipz in 1:pzNum
            pd, pz = sp.ss_pdSamples[ipd], sp.ss_pzSamples[ipz]
            spe((ti,)) = saddlePointEquation(Ax,Ay,tr,ti,pd,pz)
            ti_sol = nlsolve(spe, [asinh(ω/Ftr*sqrt(pd^2+pz^2+2Ip))/ω])
            ti = 0.0
            if converged(ti_sol) && (ti_sol.zero[1]>0)
                ti = ti_sol.zero[1]
                (ti == 0.0) && continue
            else
                continue
            end
            sample_count_thread[threadid()] += 1
            ts = tr+1im*ti   # saddle point time
            Axts, Ayts = Ax(ts), Ay(ts)
            px =  pd*imag(Ayts)/sqrt(imag(Axts)^2+imag(Ayts)^2) - real(Axts)
            py = -pd*imag(Axts)/sqrt(imag(Axts)^2+imag(Ayts)^2) - real(Ayts)
            S  = quadgk(t->((px+Ax(t))^2+(py+Ay(t))^2+pz^2)/2 + Ip, ts, tr)[1] # Integrates ∂S/∂t from ts to tr.
            x0 = quadgk(ti->imag(Ax(tr+1im*ti)), 0.0, ti)[1]
            y0 = quadgk(ti->imag(Ay(tr+1im*ti)), 0.0, ti)[1]
            z0 = 0.0
            px0 = px+Ax(tr)
            py0 = py+Ay(tr)
            pz0 = pz
            rate = exp(2*imag(S))
            init_thread[1:8,threadid(),sample_count_thread[threadid()]] = [x0,y0,z0,px0,py0,pz0,tr,rate]
            if sp.phaseMethod != :CTMC
                init_thread[9,threadid(),sample_count_thread[threadid()]] = 0.0
            end
        end
    end

    if sum(sample_count_thread) == 0
        return nothing
    end
    if sum(sample_count_thread) < pdNum*pzNum
        @warn "[SFASampleProvider] Found no root in $(pdNum*pzNum-sum(sample_count_thread)) saddle-point equations in electron batch #$batchId."
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