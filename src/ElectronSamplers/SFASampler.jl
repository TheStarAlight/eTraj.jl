
using Base.Threads
using SpecialFunctions
using StaticArrays
using Random
using NLsolve
using QuadGK

"Sample provider which generates initial electron samples through SFA formula."
struct SFASampler <: ElectronSampler
    laser   ::Laser;
    target  ::SAEAtomBase;  # SFA only supports [SAEAtomBase].
    monte_carlo;
    t_samples;
    ss_kd_samples;
    ss_kz_samples;
    mc_kt_num;
    mc_kd_max;
    mc_kz_max;
    cutoff_limit;
    phase_method;
    rate_prefix;
    mol_ion_orbit_idx;

    function SFASampler(;
        laser               ::Laser,
        target              ::SAEAtomBase,
        sample_t_intv       ::Tuple{<:Real,<:Real},
        sample_t_num        ::Integer,
        sample_cutoff_limit ::Real,
        sample_monte_carlo  ::Bool,
        traj_phase_method   ::Symbol,
        rate_prefix         ::Union{Symbol,AbstractVector{Symbol},AbstractSet{Symbol}},
            #* for step-sampling (!sample_monte_carlo)
        ss_kd_max           ::Real      = 0,
        ss_kd_num           ::Integer   = 0,
        ss_kz_max           ::Real      = 0,
        ss_kz_num           ::Integer   = 0,
            #* for Monte-Carlo-sampling (sample_monte_carlo)
        mc_kt_num           ::Integer   = 0,
        mc_kd_max           ::Real      = 0,
        mc_kz_max           ::Real      = 0,
            #* for target <: MoleculeBase
        mol_orbit_idx       ::Integer   = 0,
        kwargs...   # kwargs are surplus params.
        )

        # check phase method support.
        if ! (traj_phase_method in [:CTMC, :QTMC, :SCTS])
            error("[SFASampler] Undefined phase method [$traj_phase_method].")
            return
        end
        # check rate prefix support.
        if isa(rate_prefix, Symbol) # Exp or Full
            if rate_prefix == :Exp
                rate_prefix = []
            elseif rate_prefix == :Full
                rate_prefix = [:PreCC, :Jac]
            else
                error("[SFASampler] Undefined tunneling rate prefix [$rate_prefix].")
                return
            end
        else # a list containing Pre|PreCC, Jac.
            if length(rate_prefix) == 0
                rate_prefix = []
            elseif ! mapreduce(p->in(p,[:Pre,:PreCC,:Jac]), *, rate_prefix)
                error("[SFASampler] Undefined tunneling rate prefix [$rate_prefix].")
                return
            elseif :Pre in rate_prefix && :PreCC in rate_prefix
                error("[SFASampler] Rate prefixes [Pre] & [PreCC] conflict.")
                return
            end
        end
        # check sampling parameters.
        @assert (sample_t_num>0) "[SFASampler] Invalid time sample number $sample_t_num."
        @assert (sample_cutoff_limit≥0) "[SFASampler] Invalid cut-off limit $sample_cutoff_limit."
        if ! sample_monte_carlo # check SS sampling parameters.
            @assert (ss_kd_num>0 && ss_kz_num>0) "[SFASampler] Invalid kd/kz sample number $ss_kd_num/$ss_kz_num."
            @assert (ss_kd_max>0 && ss_kz_max>0) "[SFASampler] Invalid kd/kz sample boundaries $ss_kd_max/$ss_kz_max."
        else                    # check MC sampling parameters.
            @assert (sample_t_intv[1] < sample_t_intv[2]) "[SFASampler] Invalid sampling time interval $sample_t_intv."
            @assert (mc_kt_num>0) "[SFASampler] Invalid sampling kt_num $mc_kt_num."
            @assert (mc_kd_max>0 && mc_kz_max>0) "[ADKSampler] Invalid kd/kz sample boundaries $mc_kd_max/$mc_kz_max."
        end
        # check molecular orbital
        if typeof(target) <: MoleculeBase
            if ! (mol_orbit_idx in MolAsympCoeffAvailableIndices(target))
                MolCalcAsympCoeff!(target, mol_orbit_idx)
            end
            if MolEnergyLevels(target)[MolHOMOIndex(target)+mol_orbit_idx] ≥ 0
                error("[ADKSampler] The energy of the ionizing orbital of the molecule target is non-negative.")
            end
        end
        # finish initialization.
        return if ! sample_monte_carlo
            new(laser,target,
                sample_monte_carlo,
                range(sample_t_intv[1],sample_t_intv[2];length=sample_t_num),
                range(-abs(ss_kd_max),abs(ss_kd_max);length=ss_kd_num), range(-abs(ss_kz_max),abs(ss_kz_max);length=ss_kz_num),
                0,0,0, # for MC params. pass empty values
                sample_cutoff_limit,traj_phase_method,rate_prefix,mol_orbit_idx)
        else
            seed = 1836 # seed is mp/me :P
            t_samples = sort!(rand(MersenneTwister(seed), sample_t_num) .* (sample_t_intv[2]-sample_t_intv[1]) .+ sample_t_intv[1])
            new(laser,target,
                sample_monte_carlo,
                t_samples,
                0:0,0:0,    # for SS params. pass empty values
                mc_kt_num, mc_kd_max, mc_kz_max,
                sample_cutoff_limit,traj_phase_method,rate_prefix,mol_orbit_idx)
        end
    end
end

"Gets the total number of batches."
function batch_num(sp::SFASampler)
    return length(sp.t_samples)
end

"Gets the maximum number of electrons in a single batch. Usually the size would be smaller due to the probability cut-off."
function batch_max_size(sp::SFASampler)
    return if sp.monte_carlo
        sp.mc_kt_num
    else
        length(sp.ss_kd_samples)*length(sp.ss_kz_samples)
    end
end

"Generates a batch of electrons of `batch_id` from `sp` using SFA or MO-SFA method."
function gen_electron_batch(sp::SFASampler, batch_id::Integer)
    tr = sp.t_samples[batch_id]
    Ax::Function = LaserAx(sp.laser)
    Ay::Function = LaserAy(sp.laser)
    Fx::Function = LaserFx(sp.laser)
    Fy::Function = LaserFy(sp.laser)
    Fxtr= Fx(tr)
    Fytr= Fy(tr)
    Ftr = hypot(Fxtr, Fytr)
    ϕ   = atan(-Fytr,-Fxtr) # azimuth angle of tunneling exit
    ω   = AngFreq(sp.laser)

    target_type = if typeof(sp.target) <: SAEAtomBase
        1
    elseif typeof(sp.target) <: MoleculeBase
        2
    end
    SAEAtomTypeID = 1
    MoleculeTypeID = 2

    Z = AsympNuclCharge(sp.target)
    Ip =
        if target_type == MoleculeTypeID
            IonPotential(sp.target, sp.mol_ion_orbit_idx)
        elseif target_type == SAEAtomTypeID
            IonPotential(sp.target)
        end
    γ = AngFreq(sp.laser) * sqrt(2Ip) / (LaserF0(sp.laser) * UnitEnvelope(sp.laser)(tr)) # instantaneous Keldysh parameter
    if target_type == MoleculeTypeID
        asymp_coeff = MolAsympCoeff(sp.target, sp.mol_ion_orbit_idx)
        lMax = MolAsympCoeff_lMax(sp.target, sp.mol_ion_orbit_idx)
    elseif target_type == SAEAtomTypeID
        l = AngularQuantumNumber(sp.target)
        m = MagneticQuantumNumber(sp.target)
        C = AsympCoeff(sp.target)
    end
    prefix = sp.rate_prefix
    cutoff_limit = sp.cutoff_limit
    @inline ADKAmpExp(F,Ip,kd,kz) = exp(-(kd^2+kz^2+2*Ip)^1.5/3F)
    if Ftr == 0 || ADKAmpExp(Ftr,Ip,0.0,0.0)^2 < cutoff_limit/1e5
        return nothing
    end

    # determines Euler angles from FF to MF (α,β,γ)
    target_rot =
        if target_type == MoleculeTypeID
            MolRotation(sp.target)
        elseif target_type == SAEAtomTypeID
            θ,ϕ = QuantizationAxisOrientaion(sp.target)
            (ϕ,θ,0.0)
        end
    α,β,γ = obtain_FF_MF_Euler(target_rot, [Fxtr,Fytr])
    # computes wigner-D beforehand
    if target_type == MoleculeTypeID
        mol_wigner_D_mat = zeros(ComplexF64, lMax+1, 2lMax+1, 2lMax+1) # to get D_{m,m'}^l (α,β,γ), call wigner_D_mat[l+1, m+l+1, m_+l+1].
        for l in 0:lMax for m in -l:l for m_ in -l:l
            mol_wigner_D_mat[l+1, m+l+1, m_+l+1] = WignerD.wignerDjmn(l,m,m_, α,β,γ)
        end; end; end
    elseif target_type == SAEAtomTypeID
        atm_wigner_D_mat = zeros(ComplexF64, 2l+1, 2l+1) # to get D_{m,m'}^l (α,β,γ), call wigner_D_mat[m+l+1, m_+l+1].
        for m in -l:l for m_ in -l:l
            atm_wigner_D_mat[m+l+1, m_+l+1] = WignerD.wignerDjmn(l,m,m_, α,β,γ)
        end; end
    end
    new_x_axis, new_y_axis, new_z_axis = obtain_xyz_FF_LF(Fxtr,Fytr)

    # determines ionization amplitude (which contains phase)
    κ = sqrt(2Ip)
    n = Z/κ # n* = Z/κ
    c = 2^(n/2+1) * κ^(2n+1/2) * gamma(n/2+1) # i^((n-5)/2) is omitted, trivial
    e0 = 2.71828182845904523
    c_cc = 2^(3n/2+1) * κ^(5n+1/2) * Ftr^(-n) * (1+2γ/e0)^(-n)
    @inline px_(kd,Axts,Ayts) =  kd*imag(Ayts)/sqrt(imag(Axts)^2+imag(Ayts)^2) - real(Axts)
    @inline py_(kd,Axts,Ayts) = -kd*imag(Axts)/sqrt(imag(Axts)^2+imag(Ayts)^2) - real(Ayts)
    function saddle_point_equation(tr,ti,kd,kz)
        Axts = Ax(tr+1im*ti)
        Ayts = Ay(tr+1im*ti)
        px = px_(kd,Axts,Ayts); py = py_(kd,Axts,Ayts)
        return SVector(real(((px+Axts)^2+(py+Ayts)^2+kz^2)/2 + Ip))
    end
    function solve_spe(tr,kd,kz)
        spe((ti,)) = saddle_point_equation(tr,ti,kd,kz)
        ti_sol = nlsolve(spe, [asinh(ω/Ftr*sqrt(kd^2+kz^2+2Ip))/ω])
        ti = 0.0
        if converged(ti_sol) && (ti_sol.zero[1]>0)
            ti = ti_sol.zero[1]
        else
            return 0.0
        end
        return ti
    end
    @inline k_ts(px,py,kz,ts) = (px+Ax(ts), py+Ay(ts), kz)
    if target_type == MoleculeTypeID
        lmm_list = [(l,m,m_) for l in 0:lMax for m in -l:l for m_ in -l:l] # this cache improves the efficiency a lot
    elseif target_type == SAEAtomTypeID
        lmm_list = [(l,m,m_) for m in -l:l for m_ in -l:l]
    end
    @inline function pre_numerator((l,m,m_),kts)
        if target_type == MoleculeTypeID
            C_lm = asymp_coeff[l+1, m+l+1]
            (C_lm == 0) && return 0.0
            C_lm * mol_wigner_D_mat[l+1,m_+l+1,m+l+1] * sph_harm_lm_khat_approx(l,m_, kts..., new_x_axis, new_y_axis, new_z_axis)
        elseif target_type == SAEAtomTypeID
            C * atm_wigner_D_mat[m_+l+1,m+l+1] * sph_harm_lm_khat_approx(l,m_, kts..., new_x_axis, new_y_axis, new_z_axis)
        end
    end
    pre(kx,ky,kz,ts) = begin
        kts = k_ts(kx,ky,kz,ts)
        c * mapreduce(((l,m,m_),) -> pre_numerator((l,m,m_),kts), +, lmm_list) / (sum(kts .* (Fx(ts),Fy(ts),0.0)))^((n+1)/2)
    end
    pre_cc(kx,ky,kz,ts) = begin
        kts = k_ts(kx,ky,kz,ts)
        c_cc * mapreduce(((l,m,m_),) -> pre_numerator((l,m,m_),kts), +, lmm_list) / (sum(kts .* (Fx(ts),Fy(ts),0.0)))^((n+1)/2)
    end
    S_tun(px,py,kz,ti) = -quadgk(t->((px+Ax(t))^2+(py+Ay(t))^2+kz^2)/2+Ip, tr+1im*ti, tr)[1] # Integrates ∂S/∂t from ts to tr.
    function jac(kd,kz)
        dtr = 0.01
        dkd = 0.01
        ti_tp = solve_spe(tr+dtr,kd,kz)
        ti_tm = solve_spe(tr-dtr,kd,kz)
        ti_kp = solve_spe(tr,kd+dkd,kz)
        ti_km = solve_spe(tr,kd-dkd,kz)
        dpxdtr = px_(kd,Ax(tr+dtr+1im*ti_tp),Ay(tr+dtr+1im*ti_tp)) - px_(kd,Ax(tr-dtr+1im*ti_tm),Ay(tr-dtr+1im*ti_tm)) # actually is 4*dpx/dtr
        dpydtr = py_(kd,Ax(tr+dtr+1im*ti_tp),Ay(tr+dtr+1im*ti_tp)) - py_(kd,Ax(tr-dtr+1im*ti_tm),Ay(tr-dtr+1im*ti_tm))
        dpxdkd = px_(kd+dkd,Ax(tr+1im*ti_kp),Ay(tr+1im*ti_kp)) - px_(kd-dkd,Ax(tr+1im*ti_km),Ay(tr+1im*ti_km))
        dpydkd = py_(kd+dkd,Ax(tr+1im*ti_kp),Ay(tr+1im*ti_kp)) - py_(kd-dkd,Ax(tr+1im*ti_km),Ay(tr+1im*ti_km))
        return abs(dpxdtr*dpydkd-dpxdkd*dpydtr)/(4dtr*dkd)
    end
    step(range) = (maximum(range)-minimum(range))/length(range) # gets the step length of the range
    dkdt = if ! sp.monte_carlo
        step(sp.t_samples) * step(sp.ss_kd_samples) * step(sp.ss_kz_samples)
    else
        step(sp.t_samples) * 4*sp.mc_kd_max*sp.mc_kz_max/sp.mc_kt_num
    end
    amplitude::Function =
        if isempty(prefix)
            amp_exp(px,py,kd,kz,ti) = sqrt(dkdt) * exp(1im*S_tun(px,py,kz,ti))
        else
            if :Pre in prefix
                if :Jac in prefix
                    amp_pre_jac(px,py,kd,kz,ti) = sqrt(dkdt) * exp(1im*S_tun(px,py,kz,ti)) * pre(px,py,kz,tr+1im*ti) * sqrt(jac(kd,kz))
                else
                    amp_pre(px,py,kd,kz,ti) = sqrt(dkdt) * exp(1im*S_tun(px,py,kz,ti)) * pre(px,py,kz,tr+1im*ti)
                end
            elseif :PreCC in prefix
                if :Jac in prefix
                    amp_precc_jac(px,py,kd,kz,ti) = sqrt(dkdt) * exp(1im*S_tun(px,py,kz,ti)) * pre_cc(px,py,kz,tr+1im*ti) * sqrt(jac(kd,kz))
                else
                    amp_precc(px,py,kd,kz,ti) = sqrt(dkdt) * exp(1im*S_tun(px,py,kz,ti)) * pre_cc(px,py,kz,tr+1im*ti)
                end
            else # [:Jac]
                amp_jac(px,py,kd,kz,ti) = sqrt(dkdt) * exp(1im*S_tun(px,py,kz,ti)) * sqrt(jac(kd,kz))
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
                ti = solve_spe(tr,kd0,kz0)
                (ti == 0) && continue
                ts = tr+1im*ti
                px = px_(kd0,Ax(ts),Ay(ts))
                py = py_(kd0,Ax(ts),Ay(ts))
                x0 = quadgk(ti->imag(Ax(tr+1im*ti)), 0.0, ti)[1]
                y0 = quadgk(ti->imag(Ay(tr+1im*ti)), 0.0, ti)[1]
                z0 = 0.0
                kx0 = px+Ax(tr)
                ky0 = py+Ay(tr)
                amp = amplitude(px,py,kd0,kz0,ti)
                rate = abs2(amp)
                if rate < cutoff_limit
                    continue    # discard the sample
                end
                sample_count_thread[threadid()] += 1
                init_thread[1:8,threadid(),sample_count_thread[threadid()]] = [x0,y0,z0,kx0,ky0,kz0,tr,rate]
                if phase_method != :CTMC
                    init_thread[9,threadid(),sample_count_thread[threadid()]] = angle(amp) # Re(S_tun) is also included in the angle(amp)
                end
            end
        end
    else
        seed = 299792458
        rng = Random.MersenneTwister(seed) # use a fixed seed to ensure reproducibility
        @threads for i in 1:sp.mc_kt_num
            kd0, kz0 = gen_rand_pt(rng, sp.mc_kd_max, sp.mc_kz_max)
            ti = solve_spe(tr,kd0,kz0)
            (ti == 0) && continue
            ts = tr+1im*ti
            px = px_(kd0,Ax(ts),Ay(ts))
            py = py_(kd0,Ax(ts),Ay(ts))
            x0 = quadgk(ti->imag(Ax(tr+1im*ti)), 0.0, ti)[1]
            y0 = quadgk(ti->imag(Ay(tr+1im*ti)), 0.0, ti)[1]
            z0 = 0.0
            kx0 = px+Ax(tr)
            ky0 = py+Ay(tr)
            amp = amplitude(px,py,kd0,kz0,ti)
            rate = abs2(amp)
            if rate < cutoff_limit
                continue    # discard the sample
            end
            sample_count_thread[threadid()] += 1
            init_thread[1:8,threadid(),sample_count_thread[threadid()]] = [x0,y0,z0,kx0,ky0,kz0,tr,rate]
            if phase_method != :CTMC
                init_thread[9,threadid(),sample_count_thread[threadid()]] = angle(amp)
            end
        end
    end
    if sum(sample_count_thread) == 0
        # @warn "[SFASampler] All sampled electrons are discarded in batch #$(batchId), corresponding to t=$t."
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