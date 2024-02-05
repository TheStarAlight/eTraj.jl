
using Printf
using Base.Threads
using SpecialFunctions
using StaticArrays
using Random
using Rotations
using WignerD

"Electron sampler which generates initial electron samples using the ADK or MO-ADK method."
struct ADKSampler <: ElectronSampler
    laser   ::Laser;
    target  ::Target;
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
    adk_tun_exit;
    mol_ion_orbit_idx;

    function ADKSampler(;
        laser               ::Laser,
        target              ::Target,
        sample_t_intv       ::Tuple{<:Real,<:Real},
        sample_t_num        ::Integer,
        sample_cutoff_limit ::Real,
        sample_monte_carlo  ::Bool,
        traj_phase_method   ::Symbol,
        rate_prefix         ::Union{Symbol,AbstractVector{Symbol},AbstractSet{Symbol}},
        adk_tun_exit        ::Symbol,
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
            error("[ADKSampler] Undefined phase method [$traj_phase_method].")
            return
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
        else # a list containing Pre|PreCC, Jac.
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
        # check tunneling exit support.
        if ! (adk_tun_exit in [:IpF, :FDM, :Para])
            error("[ADKSampler] Undefined tunneling exit method [$adk_tun_exit].")
            return
        end
        # check sampling parameters.
        @assert (sample_t_num>0) "[ADKSampler] Invalid time sample number $sample_t_num."
        @assert (sample_cutoff_limit≥0) "[ADKSampler] Invalid cut-off limit $sample_cutoff_limit."
        if ! sample_monte_carlo # check SS sampling parameters.
            @assert (ss_kd_num>0 && ss_kz_num>0) "[ADKSampler] Invalid kd/kz sample number $ss_kd_num/$ss_kz_num."
            @assert (ss_kd_max>0 && ss_kz_max>0) "[ADKSampler] Invalid kd/kz sample boundaries $ss_kd_max/$ss_kz_max."
        else                    # check MC sampling parameters.
            @assert (sample_t_intv[1] < sample_t_intv[2]) "[ADKSampler] Invalid sampling time interval $sample_t_intv."
            @assert (mc_kt_num>0) "[ADKSampler] Invalid sampling kt_num $mc_kt_num."
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
        # check Keldysh parameter.
        F0 = LaserF0(laser)
        Ip = IonPotential(target)
        γ0 = AngFreq(laser) * sqrt(2Ip) / F0
        if γ0 ≥ 0.5
            @warn "[ADKSampler] Keldysh parameter γ=$(@sprintf("%.4f",γ0)), adiabatic (tunneling) condition [γ<<1] not sufficiently satisfied."
        elseif γ0 ≥ 1.0
            @warn "[ADKSampler] Keldysh parameter γ=$(@sprintf("%.4f",γ0)), adiabatic (tunneling) condition [γ<<1] unsatisfied."
        end
        # finish initialization.
        return if ! sample_monte_carlo
            new(laser,target,
                sample_monte_carlo,
                range(sample_t_intv[1],sample_t_intv[2];length=sample_t_num),
                range(-abs(ss_kd_max),abs(ss_kd_max);length=ss_kd_num), range(-abs(ss_kz_max),abs(ss_kz_max);length=ss_kz_num),
                0,0,0, # for MC params. pass empty values
                sample_cutoff_limit,traj_phase_method,rate_prefix,adk_tun_exit,mol_orbit_idx)
        else
            seed = 1836 # seed is mp/me :D
            t_samples = sort!(rand(MersenneTwister(seed), sample_t_num) .* (sample_t_intv[2]-sample_t_intv[1]) .+ sample_t_intv[1])
            new(laser,target,
                sample_monte_carlo,
                t_samples,
                0:0,0:0, # for SS params. pass empty values
                mc_kt_num, mc_kd_max, mc_kz_max,
                sample_cutoff_limit,traj_phase_method,rate_prefix,adk_tun_exit,mol_orbit_idx)
        end
    end
end

"Gets the total number of batches."
function batch_num(sp::ADKSampler)
    return length(sp.t_samples)
end

"Gets the maximum number of electrons in a single batch. Usually the size would be smaller due to the probability cut-off."
function batch_max_size(sp::ADKSampler)
    return if sp.monte_carlo
        sp.mc_kt_num
    else
        length(sp.ss_kd_samples)*length(sp.ss_kz_samples)
    end
end

"Generates a batch of electrons of `batch_id` from `sp` using ADK or MO-ADK method."
function gen_electron_batch(sp::ADKSampler, batch_id::Integer)
    tr = sp.t_samples[batch_id]
    Fxtr= LaserFx(sp.laser)(tr)
    Fytr= LaserFy(sp.laser)(tr)
    Ftr = hypot(Fxtr, Fytr)
    ϕ   = atan(-Fytr,-Fxtr) # azimuth angle of tunneling exit

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
    α,β,γ_ = obtain_FF_MF_Euler(target_rot, [Fxtr,Fytr])
    # computes wigner-D beforehand
    if target_type == MoleculeTypeID
        mol_wigner_D_mat = zeros(ComplexF64, lMax+1, 2lMax+1, 2lMax+1) # to get D_{m,m'}^l (α,β,γ), call wigner_D_mat[l+1, m+l+1, m_+l+1].
        for l in 0:lMax for m in -l:l for m_ in -l:l
            mol_wigner_D_mat[l+1, m+l+1, m_+l+1] = WignerD.wignerDjmn(l,m,m_, α,β,γ_)
        end; end; end
    elseif target_type == SAEAtomTypeID
        atm_wigner_D_mat = zeros(ComplexF64, 2l+1, 2l+1) # to get D_{m,m'}^l (α,β,γ), call wigner_D_mat[m+l+1, m_+l+1].
        for m in -l:l for m_ in -l:l
            atm_wigner_D_mat[m+l+1, m_+l+1] = WignerD.wignerDjmn(l,m,m_, α,β,γ_)
        end; end
    end
    # computes spherical harmonics beforehand
    SH_func_mat =   # to get Y_{lm}(x,y,z), call SH_func_mat[l+1,m+l+1]
        if target_type == MoleculeTypeID
            gen_sph_harm_funcs(lMax)
        elseif target_type == SAEAtomTypeID
            gen_sph_harm_funcs(l)
        end
    new_x_axis, new_y_axis, new_z_axis = obtain_xyz_FF_LF(Fxtr,Fytr)

    # determines ionization amplitude (which contains phase)
    κ = sqrt(2Ip)
    n = Z/κ # n* = Z/κ
    c = 2^(n/2+1) * κ^(2n+1/2) * gamma(n/2+1) # i^(3(1-n)/2) is omitted, trivial
    e0 = 2.71828182845904523
    c_cc = 2^(3n/2+1) * κ^(5n+1/2) * Ftr^(-n) * (1+2γ/e0)^(-n)
    @inline ti(kd,kz) = sqrt(κ^2+kd^2+kz^2)/Ftr
    @inline k_ts(kx,ky,kz,ts) = SVector{3,ComplexF64}(kx-1im*imag(ts)*Fxtr, ky-1im*imag(ts)*Fytr, kz)
    if target_type == MoleculeTypeID
        lmm_list = [(l,m,m_) for l in 0:lMax for m in -l:l for m_ in -l:l] # this cache improves the efficiency a lot
    elseif target_type == SAEAtomTypeID
        lmm_list = [(l,m,m_) for m in -l:l for m_ in -l:l]
    end
    @inline function pre_numerator((l,m,m_),kts)
        if target_type == MoleculeTypeID
            C_lm = asymp_coeff[l+1, m+l+1]
            (C_lm == 0) && return 0.0
            new_kx = kts ⋅ new_x_axis
            new_ky = kts ⋅ new_y_axis
            new_kz = kts ⋅ new_z_axis
            C_lm * mol_wigner_D_mat[l+1,m_+l+1,m+l+1] * SH_func_mat[l+1,m_+l+1](new_kx, new_ky, new_kz)
        elseif target_type == SAEAtomTypeID
            C * atm_wigner_D_mat[m_+l+1,m+l+1] * SH_func_mat[l+1,m_+l+1](new_kx, new_ky, new_kz)
        end
    end
    pre(kx,ky,kz,ts) = begin
        kts = k_ts(kx,ky,kz,ts)
        c * mapreduce(((l,m,m_),) -> pre_numerator((l,m,m_),kts), +, lmm_list) / ((kx^2+ky^2+kz^2+2Ip)*Ftr^2)^((n+1)/4)
    end
    pre_cc(kx,ky,kz,ts) = begin
        kts = k_ts(kx,ky,kz,ts)
        c_cc * mapreduce(((l,m,m_),) -> pre_numerator((l,m,m_),kts), +, lmm_list) / ((kx^2+ky^2+kz^2+2Ip)*Ftr^2)^((n+1)/4)
    end
    jac = Ftr
    step(range) = (maximum(range)-minimum(range))/length(range) # gets the step length of the range
    dkdt = if ! sp.monte_carlo
        step(sp.t_samples) * step(sp.ss_kd_samples) * step(sp.ss_kz_samples)
    else
        step(sp.t_samples) * 4*sp.mc_kd_max*sp.mc_kz_max/sp.mc_kt_num
    end
    amplitude::Function =
        if isempty(prefix)
            amp_exp(kx,ky,kd,kz,ti) = sqrt(dkdt) * ADKAmpExp(Ftr,Ip,kd,kz)
        else
            if :Pre in prefix
                if :Jac in prefix
                    amp_pre_jac(kx,ky,kd,kz,ti) = sqrt(dkdt) * ADKAmpExp(Ftr,Ip,kd,kz) * pre(kx,ky,kz,tr+1im*ti) * sqrt(jac)
                else
                    amp_pre(kx,ky,kd,kz,ti) = sqrt(dkdt) * ADKAmpExp(Ftr,Ip,kd,kz) * pre(kx,ky,kz,tr+1im*ti)
                end
            elseif :PreCC in prefix
                if :Jac in prefix
                    amp_precc_jac(kx,ky,kd,kz,ti) = sqrt(dkdt) * ADKAmpExp(Ftr,Ip,kd,kz) * pre_cc(kx,ky,kz,tr+1im*ti) * sqrt(jac)
                else
                    amp_precc(kx,ky,kd,kz,ti) = sqrt(dkdt) * ADKAmpExp(Ftr,Ip,kd,kz) * pre_cc(kx,ky,kz,tr+1im*ti)
                end
            else # [:Jac]
                amp_jac(kx,ky,kd,kz,ti) = sqrt(dkdt) * ADKAmpExp(Ftr,Ip,kd,kz) * sqrt(jac)
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
                kx0 = kd0*-sin(ϕ)
                ky0 = kd0* cos(ϕ)
                r0 = (Ip+(kd0^2+kz0^2)/2)/Ftr
                x0 = r0*cos(ϕ)
                y0 = r0*sin(ϕ)
                z0 = 0.0
                amp = amplitude(kx0,ky0,kd0,kz0,ti(kd0,kz0))
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
    else
        seed = 299792458
        rng = Random.MersenneTwister(seed) # use a fixed seed to ensure reproducibility
        @threads for i in 1:sp.mc_kt_num
            kd0, kz0 = gen_rand_pt(rng, sp.mc_kd_max, sp.mc_kz_max)
            kx0 = kd0*-sin(ϕ)
            ky0 = kd0* cos(ϕ)
            r0 = (Ip+(kd0^2+kz0^2)/2)/Ftr
            x0 = r0*cos(ϕ)
            y0 = r0*sin(ϕ)
            z0 = 0.0
            amp = amplitude(kx0,ky0,kd0,kz0,ti(kd0,kz0))
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