using StaticArrays
using NLsolve
using QuadGK
using Base.Threads

"Sample provider which yields electron samples through SFA formula, matching `IonRateMethod=:SFA`"
struct SFASampler <: ElectronSampleProvider
    laser           ::Laser;
    target          ::SAEAtomBase;          # SFA only supports [SAEAtomBase].
    tSamples        ::AbstractVector;
    ss_kdSamples    ::AbstractVector;
    ss_kzSamples    ::AbstractVector;
    phaseMethod     ::Symbol;           # currently supports :CTMC, :QTMC, :SCTS.
    ionRatePrefix   ::Symbol;           # currently supports :ExpRate.
    function SFASampler(;   laser               ::Laser,
                            target              ::SAEAtomBase,
                            sample_tSpan        ::Tuple{<:Real,<:Real},
                            sample_tSampleNum   ::Int,
                            simu_phaseMethod    ::Symbol,
                            rate_ionRatePrefix  ::Symbol,
                            ss_kdMax            ::Real,
                            ss_kdNum            ::Int,
                            ss_kzMax            ::Real,
                            ss_kzNum            ::Int,
                            kwargs...   # kwargs are surplus params.
                            )
        # check phase method support.
        if ! (simu_phaseMethod in [:CTMC, :QTMC, :SCTS])
            error("[SFASampler] Undefined phase method [$simu_phaseMethod].")
            return
        end
        # check IonRate prefix support.
        if ! (rate_ionRatePrefix in [:ExpRate, :ExpPre, :ExpJac, :Full])
            error("[SFASampler] Undefined tunneling rate prefix [$rate_ionRatePrefix].")
            return
        end
        # check sampling parameters.
        @assert (sample_tSampleNum>0) "[SFASampler] Invalid time sample number $sample_tSampleNum."
        @assert (ss_kdNum>0 && ss_kzNum>0) "[SFASampler] Invalid kd/kz sample number $ss_kdNum/$ss_kzNum."
        # finish initialization.
        return new(laser, target,
                range(sample_tSpan[1],sample_tSpan[2];length=sample_tSampleNum),
                range(-abs(ss_kdMax),abs(ss_kdMax);length=ss_kdNum), range(-abs(ss_kzMax),abs(ss_kzMax);length=ss_kzNum),
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
    Z   = AsympNuclCharge(sp.target)
    Ip  = IonPotential(sp.target)
    α   = 1.0+Z/sqrt(2Ip)
    if Ftr == 0
        return nothing
    end
    kdNum, kzNum = length(sp.ss_kdSamples), length(sp.ss_kzSamples)
    dim = (sp.phaseMethod == :CTMC) ? 8 : 9 # x,y,z,kx,ky,kz,t0,rate[,phase]
    sample_count_thread = zeros(Int,nthreads())
    init_thread = zeros(Float64, dim, nthreads(), kdNum*kzNum) # initial condition (support for multi-threading)

    "Saddle-point equation."
    function saddle_point_equation(AxF,AyF,tr,ti,kd,kz)
        Axt = AxF(tr+1im*ti)
        Ayt = AyF(tr+1im*ti)
        px =  kd*imag(Ayt)/sqrt(imag(Axt)^2+imag(Ayt)^2) - real(Axt)
        py = -kd*imag(Axt)/sqrt(imag(Axt)^2+imag(Ayt)^2) - real(Ayt)
        return SVector(real(((px+Axt)^2+(py+Ayt)^2+kz^2)/2 + Ip))
    end

    @threads for ikd in 1:kdNum
        for ikz in 1:kzNum
            kd, kz = sp.ss_kdSamples[ikd], sp.ss_kzSamples[ikz]
            spe((ti,)) = saddle_point_equation(Ax,Ay,tr,ti,kd,kz)
            ti_sol = nlsolve(spe, [asinh(ω/Ftr*sqrt(kd^2+kz^2+2Ip))/ω])
            ti = 0.0
            if converged(ti_sol) && (ti_sol.zero[1]>0)
                ti = ti_sol.zero[1]
                (ti == 0.0) && continue
            else
                continue
            end
            sample_count_thread[threadid()] += 1

            @inline function px_(AxF,AyF,tr,ti,kd)
                ts = tr+1im*ti   # saddle point time
                Axts, Ayts = AxF(ts), AyF(ts)
                return kd*imag(Ayts)/sqrt(imag(Axts)^2+imag(Ayts)^2) - real(Axts)
            end
            @inline function py_(AxF,AyF,tr,ti,kd)
                ts = tr+1im*ti   # saddle point time
                Axts, Ayts = AxF(ts), AyF(ts)
                return -kd*imag(Axts)/sqrt(imag(Axts)^2+imag(Ayts)^2) - real(Ayts)
            end
            px = px_(Ax,Ay,tr,ti,kd)
            py = py_(Ax,Ay,tr,ti,kd)
            # px =  kd*imag(Ayts)/sqrt(imag(Axts)^2+imag(Ayts)^2) - real(Axts)
            # py = -kd*imag(Axts)/sqrt(imag(Axts)^2+imag(Ayts)^2) - real(Ayts)
            S  = quadgk(t->((px+Ax(t))^2+(py+Ay(t))^2+kz^2)/2 + Ip, tr+1im*ti, tr)[1] # Integrates ∂S/∂t from ts to tr.
            x0 = quadgk(ti->imag(Ax(tr+1im*ti)), 0.0, ti)[1]
            y0 = quadgk(ti->imag(Ay(tr+1im*ti)), 0.0, ti)[1]
            z0 = 0.0
            kx0 = px+Ax(tr)
            ky0 = py+Ay(tr)
            kz0 = kz
            rate = exp(2*imag(S))
            if sp.ionRatePrefix == :ExpPre || sp.ionRatePrefix == :Full
                rate /= abs((px+Ax(tr+1im*ti))*Fx(tr+1im*ti)+(py+Ay(tr+1im*ti))*Fy(tr+1im*ti))^α
            end
            if sp.ionRatePrefix == :ExpJac || sp.ionRatePrefix == :Full
                dt = 0.01
                dk = 0.01
                spe_t_p_dt((ti,)) = saddle_point_equation(Ax,Ay,tr+dt,ti,kd,kz)
                spe_t_m_dt((ti,)) = saddle_point_equation(Ax,Ay,tr-dt,ti,kd,kz)
                spe_k_p_dk((ti,)) = saddle_point_equation(Ax,Ay,tr,ti,kd+dk,kz)
                spe_k_m_dk((ti,)) = saddle_point_equation(Ax,Ay,tr,ti,kd-dk,kz)
                ti_tp = nlsolve(spe_t_p_dt,[ti]).zero[1]
                ti_tm = nlsolve(spe_t_m_dt,[ti]).zero[1]
                ti_kp = nlsolve(spe_k_p_dk,[ti]).zero[1]
                ti_km = nlsolve(spe_k_m_dk,[ti]).zero[1]
                dpxdtr = px_(Ax,Ay,tr+dt,ti_tp,kd) - px_(Ax,Ay,tr-dt,ti_tm,kd)
                dpydtr = py_(Ax,Ay,tr+dt,ti_tp,kd) - py_(Ax,Ay,tr-dt,ti_tm,kd)
                dpxdkd = px_(Ax,Ay,tr,ti_kp,kd+dk) - px_(Ax,Ay,tr,ti_km,kd-dk)
                dpydkd = py_(Ax,Ay,tr,ti_kp,kd+dk) - py_(Ax,Ay,tr,ti_km,kd-dk)
                rate *= abs(dpxdtr*dpydkd-dpxdkd*dpydtr)/(4dt*dk)
            end
            init_thread[1:8,threadid(),sample_count_thread[threadid()]] = [x0,y0,z0,kx0,ky0,kz0,tr,rate]
            if sp.phaseMethod != :CTMC
                init_thread[9,threadid(),sample_count_thread[threadid()]] = 0.0
            end
        end
    end

    if sum(sample_count_thread) == 0
        return nothing
    end
    if sum(sample_count_thread) < kdNum*kzNum
        @warn "[SFASampler] Found no root in $(kdNum*kzNum-sum(sample_count_thread)) saddle-point equations in electron batch #$batchId."
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