
using Printf
using Base.Threads
using Random
using Rotations
using WignerD
using HDF5

"Electron sampler which generates initial electron samples using the WFAT formula."
struct WFATSampler <: ElectronSampler
    laser   ::Laser;
    target  ::MoleculeBase; # WFAT only supports [MoleculeBase]
    monte_carlo;
    t_samples;
    ss_kd_samples;
    ss_kz_samples;
    mc_kt_num;
    mc_kd_max;
    mc_kz_max;
    cutoff_limit;
    phase_method;   # currently supports :CTMC.
    mol_ion_orbit_idx;

    function WFATSampler(;
        laser               ::Laser,
        target              ::MoleculeBase,
        sample_t_intv       ::Tuple{<:Real,<:Real},
        sample_t_num        ::Integer,
        sample_cutoff_limit ::Real,
        sample_monte_carlo  ::Bool,
        traj_phase_method   ::Symbol,
        mol_orbit_idx       ::Integer,
            #* for step-sampling (!sample_monte_carlo)
        ss_kd_max           ::Real      = 0,
        ss_kd_num           ::Integer   = 0,
        ss_kz_max           ::Real      = 0,
        ss_kz_num           ::Integer   = 0,
            #* for Monte-Carlo-sampling (sample_monte_carlo)
        mc_kt_num           ::Integer   = 0,
        mc_kd_max           ::Real      = 0,
        mc_kz_max           ::Real      = 0,
        kwargs...   # kwargs are surplus params.
        )

        # check phase method support.
        if ! (traj_phase_method in [:CTMC])
            error("[WFATSampler] Unsupported phase method [$traj_phase_method].")
            return
        end
        # check sampling parameters.
        @assert (sample_t_num>0) "[WFATSampler] Invalid time sample number $sample_t_num."
        @assert (sample_cutoff_limit≥0) "[WFATSampler] Invalid cut-off limit $sample_cutoff_limit."
        if ! sample_monte_carlo # check SS sampling parameters.
            @assert (ss_kd_num>0 && ss_kz_num>0) "[WFATSampler] Invalid kd/kz sample number $ss_kd_num/$ss_kz_num."
        else                    # check MC sampling parameters.
            @assert (sample_t_intv[1] < sample_t_intv[2]) "[WFATSampler] Invalid sampling time interval $sample_t_intv."
            @assert (mc_kt_num>0) "[WFATSampler] Invalid sampling kt_num $mc_kt_num."
            @assert (mc_kd_max>0 && mc_kz_max>0) "[WFATSampler] Invalid kd/kz sample boundaries $mc_kd_max/$mc_kz_max."
        end
        # check molecular orbital data
        if ! (mol_orbit_idx in MolWFATAvailableIndices(target))
            MolCalcWFATData!(target, mol_orbit_idx)
        end
        if MolEnergyLevels(target)[MolHOMOIndex(target)+mol_orbit_idx] ≥ 0
            error("[WFATSampler] The energy of the ionizing orbit is non-negative.")
        end
        # check Keldysh parameter.
        if laser isa MonochromaticLaser
            γ0 = KeldyshParameter(laser, IonPotential(target))
            if γ0 ≥ 0.5
                @warn "[WFATSampler] Keldysh parameter γ=$(@sprintf("%.4f",γ0)), adiabatic (tunneling) condition [γ<<1] not sufficiently satisfied."
            elseif γ0 ≥ 1.0
                @warn "[WFATSampler] Keldysh parameter γ=$(@sprintf("%.4f",γ0)), adiabatic (tunneling) condition [γ<<1] unsatisfied."
            end
        elseif laser isa BichromaticLaser
            γ10 = KeldyshParameter(laser[1], IonPotential(target))
            γ20 = KeldyshParameter(laser[2], IonPotential(target))
            if max(γ10,γ20) ≥ 0.5
                @warn "[WFATSampler] Keldysh parameter γ₁=$(@sprintf("%.4f",γ10)), γ₂=$(@sprintf("%.4f",γ20)), adiabatic (tunneling) condition [γ<<1] not sufficiently satisfied."
            elseif max(γ10,γ20) ≥ 1.0
                @warn "[WFATSampler] Keldysh parameter γ₁=$(@sprintf("%.4f",γ10)), γ₂=$(@sprintf("%.4f",γ20)), adiabatic (tunneling) condition [γ<<1] unsatisfied."
            end
        end
        # finish initialization.
        return if ! sample_monte_carlo
            new(laser,target,
                sample_monte_carlo,
                range(sample_t_intv[1],sample_t_intv[2];length=sample_t_num),
                range(-abs(ss_kd_max),abs(ss_kd_max);length=ss_kd_num), range(-abs(ss_kz_max),abs(ss_kz_max);length=ss_kz_num),
                0,0,0,  # for MC params. pass empty values
                sample_cutoff_limit,traj_phase_method,mol_orbit_idx)
        else
            seed = 1836 # seed is mp/me
            t_samples = sort!(rand(MersenneTwister(seed), sample_t_num) .* (sample_t_intv[2]-sample_t_intv[1]) .+ sample_t_intv[1])
            new(laser,target,
                sample_monte_carlo,
                t_samples,
                0:0,0:0,    # for SS params. pass empty values
                mc_kt_num, mc_kd_max, mc_kz_max,
                sample_cutoff_limit,traj_phase_method,mol_orbit_idx)
        end
    end
end

"Gets the total number of batches."
function batch_num(sp::WFATSampler)
    return length(sp.t_samples)
end

"Gets the maximum number of electrons in a single batch. Usually the size would be smaller due to the probability cut-off."
function batch_max_size(sp::WFATSampler)
    return if sp.monte_carlo
        sp.mc_kt_num
    else
        length(sp.ss_kd_samples)*length(sp.ss_kz_samples)
    end
end

"Generates a batch of electrons of `batch_id` from `sp` using WFAT method."
function gen_electron_batch(sp::WFATSampler, batch_id::Integer)
    tr = sp.t_samples[batch_id]
    Fxtr = LaserFx(sp.laser)(tr)
    Fytr = LaserFy(sp.laser)(tr)
    Ftr = hypot(Fxtr,Fytr)
    ϕ  = atan(-Fytr,-Fxtr)   # direction of tunneling exit, which is opposite to F.
    Z  = AsympNuclCharge(sp.target)
    Ip = IonPotential(sp.target, sp.mol_ion_orbit_idx)
    @inline ADKAmpExp(F,Ip,kd,kz) = exp(-(kd^2+kz^2+2*Ip)^1.5/3F)
    cutoff_limit = sp.cutoff_limit
    if Ftr == 0 || ADKAmpExp(Ftr,Ip,0.0,0.0)^2 < cutoff_limit/1e5
        return nothing
    end

    # determines Euler angles from FF to MF (α,β,γ)
    mol_rot = MolRotation(sp.target)
    α,β,γ = obtain_FF_MF_Euler(mol_rot, [Fxtr,Fytr])

    # determines ionization rate
    κ = sqrt(2Ip)
    n = Z/κ # n* = Z/κ
    # prepare G (structure factor)
    nξMax, mMax = MolWFATMaxChannels(sp.target, sp.mol_ion_orbit_idx)
    G_data = zeros(ComplexF64, nξMax+1, 2mMax+1)    # to obtain G_nξ,m, call G_data[nξ+1, m+mMax+1]
    for nξ in 0:nξMax, m in -mMax:mMax
        G_data[nξ+1, m+mMax+1] = MolWFATStructureFactor_G(sp.target,sp.mol_ion_orbit_idx,nξ,m,β,γ)
    end
    # define W_F (modified field factor)
    W(nξ,abs_m,kd,kz) = (κ^(abs_m+2)/2/Ftr^(abs_m+1)/factorial(abs_m)) * (4κ^2/Ftr)^(2n-2nξ-abs_m-1) * (kd^2+kz^2)^abs_m * exp(-2*(κ^2+kd^2+kz^2)^1.5/3Ftr)
    # calc dkdt
    step(range) = (maximum(range)-minimum(range))/length(range) # gets the step length of the range
    dkdt = if ! sp.monte_carlo
        step(sp.t_samples) * step(sp.ss_kd_samples) * step(sp.ss_kz_samples)
    else
        step(sp.t_samples) * 4*sp.mc_kd_max*sp.mc_kz_max/sp.mc_kt_num
    end
    rate(kd,kz) = sum(((nξ,m),)->abs2(G_data[nξ+1,m+mMax+1])*W(nξ,abs(m),kd,kz), [(nξ,m) for nξ in 0:nξMax, m in -mMax:mMax]) * dkdt

    dim = 8 # x,y,z,kx,ky,kz,t0,rate
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
                rate_ = rate(kd0,kz0)
                if rate_ < cutoff_limit
                    continue    # discard the sample
                end
                sample_count_thread[threadid()] += 1
                init_thread[1:8,threadid(),sample_count_thread[threadid()]] = [x0,y0,z0,kx0,ky0,kz0,tr,rate_]
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
            rate_ = rate(kd0,kz0)
            if rate_ < cutoff_limit
                continue    # discard the sample
            end
            sample_count_thread[threadid()] += 1
            init_thread[1:8,threadid(),sample_count_thread[threadid()]] = [x0,y0,z0,kx0,ky0,kz0,tr,rate_]
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

