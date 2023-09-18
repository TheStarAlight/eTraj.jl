
using Printf
using Base.Threads
using SpecialFunctions
using StaticArrays
using Random
"Sample provider which generates initial electron samples using the ADK method."
struct ADKSampler <: ElectronSampleProvider
    laser   ::Laser;
    target  ::SAEAtomBase;  # ADK only supports [SAEAtomBase].
    monte_carlo;
    t_samples;
    ss_kd_samples;
    ss_kz_samples;
    mc_kt_num;
    mc_kt_max;
    cutoff_limit;
    phase_method;   # currently supports :CTMC, :QTMC, :SCTS.
    rate_prefix;    # supports :Exp, :Full or a combination of {:Pre|:PreCC, :Jac}.
    ADK_tun_exit;   # currently supports :IpF, :FDM, :Para.

    function ADKSampler(;
                        laser               ::Laser,
                        target              ::SAEAtomBase,
                        sample_t_intv       ::Tuple{<:Real,<:Real},
                        sample_t_num        ::Int,
                        sample_cutoff_limit ::Real,
                        sample_monte_carlo  ::Bool,
                        traj_phase_method   ::Symbol,
                        rate_prefix         ::Union{Symbol,AbstractVector{Symbol},AbstractSet{Symbol}},
                        adk_tun_exit        ::Symbol,
                            #* for step-sampling (!sample_monte_carlo)
                        ss_kd_max           ::Real,
                        ss_kd_num           ::Int,
                        ss_kz_max           ::Real,
                        ss_kz_num           ::Int,
                            #* for Monte-Carlo-sampling (sample_monte_carlo)
                        mc_kt_num           ::Int,
                        mc_kt_max           ::Real,
                        kwargs...   # kwargs are surplus params.
                        )
        F0 = LaserF0(laser)
        Ip = IonPotential(target)
        γ0 = AngFreq(laser) * sqrt(2Ip) / F0
        # check phase method support.
        if ! (traj_phase_method in [:CTMC, :QTMC, :SCTS])
            error("[ADKSampler] Undefined phase method [$traj_phase_method].")
            return
        end
        # check tunneling exit support.
        if ! (adk_tun_exit in [:IpF, :FDM, :Para])
            error("[ADKSampler] Undefined tunneling exit method [$adk_tun_exit].")
            return
        end
        # switch tunneling exit method if inappropriate.
        if adk_tun_exit == :FDM && Ip^2 < 4F0
            adk_tun_exit == :IpF
            @warn "[ADKSampler] Tunneling exit method has changed from [FDM] to [IpF] due to failure in tunneling exit determination of [FDM] model."
        elseif adk_tun_exit == :Para && Ip^2 < 4*(1-sqrt(Ip/2))*F0
            adk_tun_exit == :IpF
            @warn "[ADKSampler] Tunneling exit method has changed from [Para] to [IpF] due to failure in tunneling exit determination of [Para] model."
        end
        # check rate prefix support.
        if isa(rate_prefix, Symbol) # Exp or Full
            if rate_prefix == :Exp
                rate_prefix = []
            elseif rate_prefix == :Full
                rate_prefix = [:PreCC, :Jac]
            else
                error("[ADKSampler] Undefined tunneling rate prefix [$rate_prefix].")
                return
            end
        else    # a list containing Pre|PreCC, Jac.
            if length(rate_prefix) == 0
                rate_prefix = []
            elseif ! mapreduce(p->in(p,[:Pre,:PreCC,:Jac]), *, rate_prefix)
                error("[ADKSampler] Undefined tunneling rate prefix [$rate_prefix].")
                return
            elseif :Pre in rate_prefix && :PreCC in rate_prefix
                error("[ADKSampler] Rate prefixes [Pre] & [PreCC] conflict.")
                return
            end
        end
        # check Keldysh parameter.
        if γ0 ≥ 0.5
            @warn "[ADKSampler] Keldysh parameter γ=$(@sprintf "%.4f" γ0), adiabatic (tunneling) condition [γ<<1] not sufficiently satisfied."
        elseif γ0 ≥ 1.0
            @warn "[ADKSampler] Keldysh parameter γ=$(@sprintf "%.4f" γ0), adiabatic (tunneling) condition [γ<<1] unsatisfied."
        end
        # check sampling parameters.
        @assert (sample_t_num>0) "[ADKSampler] Invalid time sample number $sample_t_num."
        @assert (sample_cutoff_limit≥0) "[ADKSampler] Invalid cut-off limit $sample_cutoff_limit."
        if ! sample_monte_carlo # check SS sampling parameters.
            @assert (ss_kd_num>0 && ss_kz_num>0) "[ADKSampler] Invalid kd/kz sample number $ss_kd_num/$ss_kz_num."
        else                    # check MC sampling parameters.
            @assert (sample_t_intv[1] < sample_t_intv[2]) "[ADKSampler] Invalid sampling time interval $sample_t_intv."
            @assert (mc_kt_num>0) "[ADKSampler] Invalid sampling kt_num $mc_kt_num."
            @assert (mc_kt_max>0) "[ADKSampler] Invalid sampling kt_max $mc_kt_max."
        end
        # finish initialization.
        return if ! sample_monte_carlo
            new(laser,target,
                sample_monte_carlo,
                range(sample_t_intv[1],sample_t_intv[2];length=sample_t_num),
                range(-abs(ss_kd_max),abs(ss_kd_max);length=ss_kd_num), range(-abs(ss_kz_max),abs(ss_kz_max);length=ss_kz_num),
                0,0,        # for MC params. pass empty values
                sample_cutoff_limit,traj_phase_method,rate_prefix,adk_tun_exit)
        else
            t_samples = sort!(rand(MersenneTwister(1), sample_t_num) .* (sample_t_intv[2]-sample_t_intv[1]) .+ sample_t_intv[1])
            new(laser,target,
                sample_monte_carlo,
                t_samples,
                0:0,0:0,    # for SS params. pass empty values
                mc_kt_num, mc_kt_max,
                sample_cutoff_limit,traj_phase_method,rate_prefix,adk_tun_exit)
        end
    end
end
"Gets the total number of batches."
function batch_num(sp::ADKSampler)
    return length(sp.t_samples)
end

"Generates a batch of electrons of `batchId` from `sp` using ADK method."
function gen_electron_batch(sp::ADKSampler, batchId::Integer)
    t = sp.t_samples[batchId]
    Fx::Function = LaserFx(sp.laser)
    Fy::Function = LaserFy(sp.laser)
    Fxt = Fx(t)
    Fyt = Fy(t)
    Ft  = hypot(Fxt,Fyt)
    F0  = LaserF0(sp.laser)
    φ   = atan(-Fyt,-Fxt)
    Z   = AsympNuclCharge(sp.target)
    Ip  = IonPotential(sp.target)
    l   = AngularQuantumNumber(sp.target)
    m   = MagneticQuantumNumber(sp.target)
    C   = AsympCoeff(sp.target)
    γ0  = AngFreq(sp.laser) * sqrt(2Ip) / F0
    prefix = sp.rate_prefix
    @inline ADKAmpExp(F,Ip,kd,kz) = exp(-(kd^2+kz^2+2*Ip)^1.5/3F)
    cutoff_limit = sp.cutoff_limit
    if Ft == 0 || ADKAmpExp(Ft,Ip,0.0,0.0)^2 < cutoff_limit/1e3
        return nothing
    end
    # determining tunneling exit position. here Ip_eff=Ip+k0^2/2
    r_exit =
        if      sp.ADK_tun_exit == :IpF
            Ip_eff -> Ip_eff / Ft
        elseif  sp.ADK_tun_exit == :FDM
            Ip_eff -> (Ip_eff + sqrt(Ip_eff^2 - 4Ft*Z)) / 2Ft
        elseif  sp.ADK_tun_exit == :Para
            Ip_eff -> (Ip_eff + sqrt(Ip_eff^2 - 4*(Z-sqrt(Ip_eff/2)*(1+abs(m)))*Ft)) / 2Ft
        end
    # determining ionization amplitude (contains phase)
    amplitude::Function =
        begin
            κ = sqrt(2Ip)
            n = Z/κ # n* = Z/κ
            c = 2^(n/2+1) * κ^(2n+1/2) * gamma(n/2+1) # i^(3(1-n)/2) is omitted, trivial
            e0 = 2.71828182845904523
            c_cc = 2^(3n/2+1) * κ^(5n+1/2) * Ft^(-n) * (1+2γ0/e0)^(-n)
            @inline ti(kd,kz) = sqrt(κ^2+kd^2+kz^2)/Ft
            @inline k_ts(kx,ky,kz,ts) = (kx,ky,kz) .- (1im*imag(ts) .* (Fxt,Fyt,0.0))
            new_x_axis = @SVector [0.0,0.0,1.0]
            new_z_axis = @SVector [Fxt/Ft,Fyt/Ft,0.0]
            new_y_axis = new_z_axis × new_x_axis
            pre(kx,ky,kz,ts) = begin
                kts = k_ts(kx,ky,kz,ts)
                c * C * sph_harm_lm_khat(l,m, kts..., new_x_axis, new_y_axis, new_z_axis) / ((kx^2+ky^2+kz^2+2Ip)*Ft^2)^((n+1)/4)
            end
            pre_cc(kx,ky,kz,ts) = begin
                kts = k_ts(kx,ky,kz,ts)
                c_cc * C * sph_harm_lm_khat(l,m, kts..., new_x_axis, new_y_axis, new_z_axis) / ((kx^2+ky^2+kz^2+2Ip)*Ft^2)^((n+1)/4)
            end
            jac = Ft
            step(range) = (maximum(range)-minimum(range))/length(range) # gets the step length of the range
            dkdt = if ! sp.monte_carlo
                step(sp.t_samples) * step(sp.ss_kd_samples) * step(sp.ss_kz_samples)
            else
                step(sp.t_samples) * π*sp.mc_kt_max^2/sp.mc_kt_num
            end

            # returns
            if isempty(prefix)
                amp_exp(kx,ky,kd,kz,ti) = sqrt(dkdt) * ADKAmpExp(Ft,Ip,kd,kz)
            else
                if :Pre in prefix
                    if :Jac in prefix
                        amp_pre_jac(kx,ky,kd,kz,ti) = sqrt(dkdt) * ADKAmpExp(Ft,Ip,kd,kz) * pre(kx,ky,kz,t+1im*ti) * sqrt(jac)
                    else
                        amp_pre(kx,ky,kd,kz,ti) = sqrt(dkdt) * ADKAmpExp(Ft,Ip,kd,kz) * pre(kx,ky,kz,t+1im*ti)
                    end
                elseif :PreCC in prefix
                    if :Jac in prefix
                        amp_precc_jac(kx,ky,kd,kz,ti) = sqrt(dkdt) * ADKAmpExp(Ft,Ip,kd,kz) * pre_cc(kx,ky,kz,t+1im*ti) * sqrt(jac)
                    else
                        amp_precc(kx,ky,kd,kz,ti) = sqrt(dkdt) * ADKAmpExp(Ft,Ip,kd,kz) * pre_cc(kx,ky,kz,t+1im*ti)
                    end
                else # [:Jac]
                    amp_jac(kx,ky,kd,kz,ti) = sqrt(dkdt) * ADKAmpExp(Ft,Ip,kd,kz) * sqrt(jac)
                end
            end
        end
    phase_method = sp.phase_method
    dim = (phase_method == :CTMC) ? 8 : 9 # x,y,z,kx,ky,kz,t0,rate[,phase]

    sample_count_thread = zeros(Int,nthreads())
    init_thread = if ! sp.monte_carlo
        zeros(Float64, dim, nthreads(), length(sp.ss_kd_samples)*length(sp.ss_kz_samples)) # initial condition (support for multi-threading)
    else
        zeros(Float64, dim, nthreads(), sp.mc_kt_num)
    end

    if ! sp.monte_carlo
        kd_samples = sp.ss_kd_samples
        kz_samples = sp.ss_kz_samples
        kdNum, kzNum = length(kd_samples), length(kz_samples)
        @threads for ikd in 1:kdNum
            for ikz in 1:kzNum
                kd0, kz0 = kd_samples[ikd], kz_samples[ikz]
                kx0 = kd0*-sin(φ)
                ky0 = kd0* cos(φ)
                r0 = r_exit(Ip+(kd0^2+kz0^2)/2)
                x0 = r0*cos(φ)
                y0 = r0*sin(φ)
                z0 = 0.0
                amp = amplitude(kx0,ky0,kd0,kz0,ti(kd0,kz0))
                rate = abs2(amp)
                if rate < cutoff_limit
                    continue    # discard the sample
                end
                sample_count_thread[threadid()] += 1
                init_thread[1:8,threadid(),sample_count_thread[threadid()]] = [x0,y0,z0,kx0,ky0,kz0,t,rate]
                if phase_method != :CTMC
                    init_thread[9,threadid(),sample_count_thread[threadid()]] = angle(amp)
                end
            end
        end
    else
        rng = Random.MersenneTwister(0) # use a fixed seed to ensure reproducibility
        @threads for i in 1:sp.mc_kt_num
            # generates random (kd0,kz0) inside circle kd0^2+kz0^2=ktMax^2.
            kd0, kz0 = gen_rand_pt_circ(rng, sp.mc_kt_max)
            kx0 = kd0*-sin(φ)
            ky0 = kd0* cos(φ)
            r0 = r_exit(Ip+(kd0^2+kz0^2)/2)
            x0 = r0*cos(φ)
            y0 = r0*sin(φ)
            z0 = 0.0
            amp = amplitude(kx0,ky0,kd0,kz0,ti(kd0,kz0))
            rate = abs2(amp)
            if rate < cutoff_limit
                continue    # discard the sample
            end
            sample_count_thread[threadid()] += 1
            init_thread[1:8,threadid(),sample_count_thread[threadid()]] = [x0,y0,z0,kx0,ky0,kz0,t,rate]
            if phase_method != :CTMC
                init_thread[9,threadid(),sample_count_thread[threadid()]] = angle(amp)
            end
        end
    end
    if sum(sample_count_thread) == 0
        # @warn "[ADKSampler] All sampled electrons are discarded in batch #$(batchId), corresponding to t=$t."
        return nothing
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