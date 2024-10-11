
"Electron sampler which generates initial electron samples using the SFA-SPANE or MO-SFA-SPANE method."
struct SPANESampler <: ElectronSampler
    laser   ::Laser;
    target  ::Union{SAEAtomBase, MoleculeBase};
    dimension;
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
    mol_orbit_ridx;

    function SPANESampler(;
        laser               ::Laser,
        target              ::Union{SAEAtomBase, MoleculeBase},
        dimension           ::Integer,
        sample_t_intv       ::Tuple{<:Real,<:Real},
        sample_t_num        ::Integer,
        sample_cutoff_limit ::Real,
        sample_monte_carlo  ::Bool,
        traj_phase_method   ::Symbol,
        rate_prefix         ::Union{Symbol,AbstractVector{Symbol},AbstractSet{Symbol}},
            #* for step-sampling (!sample_monte_carlo)
        ss_kd_max           ::Real,
        ss_kd_num           ::Integer,
        ss_kz_max           ::Real,
        ss_kz_num           ::Integer,
            #* for Monte-Carlo-sampling (sample_monte_carlo)
        mc_kt_num           ::Integer,
        mc_kd_max           ::Real,
        mc_kz_max           ::Real,
            #* for target <: MoleculeBase
        mol_orbit_ridx,
        kwargs...   # kwargs are surplus params.
        )

        # check phase method support.
        if ! (traj_phase_method in [:CTMC, :QTMC, :SCTS])
            error("[SPANESampler] Undefined phase method [$traj_phase_method].")
            return
        end
        # check rate prefix support.
        if rate_prefix isa Symbol # Exp or Full
            if rate_prefix == :Exp
                rate_prefix = Set()
            elseif rate_prefix == :Full
                rate_prefix = Set([:PreCC, :Jac])
            elseif rate_prefix in [:Pre, :PreCC, :Jac]
                rate_prefix = Set([rate_prefix])
            else
                error("[SPANESampler] Undefined tunneling rate prefix [$rate_prefix].")
                return
            end
        else # a collection containing Pre|PreCC, Jac.
            if length(rate_prefix) == 0
                rate_prefix = []
            elseif ! mapreduce(p->in(p,[:Pre,:PreCC,:Jac]), *, rate_prefix)
                error("[SPANESampler] Undefined tunneling rate prefix [$rate_prefix].")
                return
            elseif :Pre in rate_prefix && :PreCC in rate_prefix
                error("[SPANESampler] Rate prefixes [Pre] & [PreCC] conflict.")
                return
            end
            rate_prefix = Set(rate_prefix)
        end
        if :PreCC in rate_prefix && !(laser isa MonochromaticLaser)
            @warn "[SPANESampler] Laser is not monochromatic, Coulomb correction in rate prefix is unavailable."
            # replace :PreCC with :Pre
            delete!(rate_prefix, :PreCC)
            push!(rate_prefix, :Pre)
        end
        # check molecular orbital
        if target isa MoleculeBase
            if ! (mol_orbit_ridx in MolAsympCoeffAvailableIndices(target))
                MolCalcAsympCoeff!(target, mol_orbit_ridx)
            end
            if IonPotential(target, mol_orbit_ridx) <= 0
                error("[SPANESampler] The energy of the ionizing orbital of the molecule target is non-negative.")
            end
        end
        # check Keldysh paramater.
        if laser isa MonochromaticLaser
            γ0 = KeldyshParameter(laser, IonPotential(target))
            if γ0 ≥ 1.0
                @warn "[SPANESampler] Keldysh parameter γ=$(@sprintf("%.4f",γ0)), adiabatic (tunneling) condition [γ<<1] unsatisfied."
            end
        elseif laser isa BichromaticLaser
            γ10 = KeldyshParameter(laser[1], IonPotential(target))
            γ20 = KeldyshParameter(laser[2], IonPotential(target))
            if max(γ10,γ20) ≥ 1.0
                @warn "[SPANESampler] Keldysh parameter γ₁=$(@sprintf("%.4f",γ10)), γ₂=$(@sprintf("%.4f",γ20)), adiabatic (tunneling) condition [γ<<1] unsatisfied."
            end
        end
        # check sampling parameters.
        @assert (sample_t_num>0) "[SPANESampler] Invalid time sample number $sample_t_num."
        @assert (sample_cutoff_limit≥0) "[SPANESampler] Invalid cut-off limit $sample_cutoff_limit."
        if dimension == 3
            if ! sample_monte_carlo
                @assert (ss_kd_num>0 && ss_kz_num>0) "[SPANESampler] Invalid kd,kz sample number $ss_kd_num,$ss_kz_num."
                @assert (ss_kd_max>0 || ss_kz_max>0) "[SPANESampler] Invalid kd,kz sample boundaries $ss_kd_max,$ss_kz_max (cannot be zero simultaneously)."
                kd_samples = range(-abs(ss_kd_max),abs(ss_kd_max); length=ss_kd_num)
                kz_samples = range(-abs(ss_kz_max),abs(ss_kz_max); length=ss_kz_num)
            else
                @assert (sample_t_intv[1] < sample_t_intv[2]) "[SPANESampler] Invalid sampling time interval $sample_t_intv."
                @assert mc_kt_num>0 "[SPANESampler] Invalid sampling kt_num $mc_kt_num."
                @assert (mc_kd_max>0 || mc_kz_max>0) "[SPANESampler] Invalid kd,kz sample boundaries $mc_kd_max,$mc_kz_max (cannot be zero simultaneously)."
            end
        else # dimension == 2
            if ! sample_monte_carlo
                @assert ss_kd_num>0 "[SPANESampler] Invalid kd sample number $ss_kd_num."
                @assert ss_kd_max>0 "[SPANESampler] Invalid kd sample boundaries $ss_kd_max."
                ss_kz_num, ss_kz_max = 1, 0.0
                kd_samples = range(-abs(ss_kd_max),abs(ss_kd_max); length=ss_kd_num)
                kz_samples = range(-abs(ss_kz_max),abs(ss_kz_max); length=ss_kz_num)
            else
                @assert (sample_t_intv[1] < sample_t_intv[2]) "[SPANESampler] Invalid sampling time interval $sample_t_intv."
                @assert mc_kt_num>0 "[SPANESampler] Invalid sampling kt_num $mc_kt_num."
                @assert mc_kd_max>0 "[SPANESampler] Invalid kd sample boundaries $mc_kd_max."
                ss_kz_max = 0.0
            end
        end
        # finish initialization.
        return if ! sample_monte_carlo
            new(laser,target,dimension,
                sample_monte_carlo,
                range(sample_t_intv[1],sample_t_intv[2];length=sample_t_num),
                kd_samples, kz_samples,
                0,0,0, # for MC params. pass empty values
                sample_cutoff_limit,traj_phase_method,rate_prefix,mol_orbit_ridx)
        else
            t_samples = sort!(rand(sample_t_num) .* (sample_t_intv[2]-sample_t_intv[1]) .+ sample_t_intv[1])
            new(laser,target,dimension,
                sample_monte_carlo,
                t_samples,
                0:0,0:0,    # for SS params. pass empty values
                mc_kt_num, mc_kd_max, mc_kz_max,
                sample_cutoff_limit,traj_phase_method,rate_prefix,mol_orbit_ridx)
        end
    end
end

"Gets the total number of batches."
function batch_num(sp::SPANESampler)
    return length(sp.t_samples)
end

"Gets the maximum number of electrons in a single batch. Usually the size would be smaller due to the probability cut-off."
function batch_max_size(sp::SPANESampler)
    return if sp.monte_carlo
        sp.mc_kt_num
    else
        length(sp.ss_kd_samples)*length(sp.ss_kz_samples)
    end
end

"Generates a batch of electrons of `batch_id` from `sp` using SFA-SPANE or MO-SFA-SPANE method."
function gen_electron_batch(sp::SPANESampler, batch_id::Integer)
    target = sp.target
    tr = sp.t_samples[batch_id]
    Fx::Function = LaserFx(sp.laser)
    Fy::Function = LaserFy(sp.laser)
    dFxt = derivative(Fx,tr)
    dFyt = derivative(Fy,tr)
    Fxtr = Fx(tr)
    Fytr = Fy(tr)
    Ftr  = hypot(Fxtr,Fytr)
    φ    = atan(-Fytr,-Fxtr)

    Z = AsympNuclCharge(sp.target)
    Ip =
        if target isa MoleculeBase
            IonPotential(sp.target, sp.mol_orbit_ridx)
        elseif target isa SAEAtomBase
            IonPotential(sp.target)
        end
    γ = if sp.laser isa MonochromaticLaser
        AngFreq(sp.laser) * sqrt(2Ip) / (LaserF0(sp.laser) * UnitEnvelope(sp.laser)(tr)) # instantaneous Keldysh parameter
    else    # cannot define γ when laser is not monochromatic
        0.0
    end
    if target isa MoleculeBase
        asymp_coeff = MolAsympCoeff(sp.target, sp.mol_orbit_ridx)
        lMax = MolAsympCoeff_lMax(sp.target, sp.mol_orbit_ridx)
    elseif target isa SAEAtomBase
        l = AngularQuantumNumber(sp.target)
        m = MagneticQuantumNumber(sp.target)
        C = AsympCoeff(sp.target)
    end
    prefix = sp.rate_prefix
    cutoff_limit = sp.cutoff_limit
    @inline ADKAmpExp(F,Ip,kd,kz) = exp(-(kd^2+kz^2+2*Ip)^1.5/3F)

    # determines Euler angles from FF to MF (α,β,γ)
    target_rot =
        if target isa MoleculeBase
            MolRotation(target)
        elseif target isa SAEAtomBase
            θ,ϕ = QuantizationAxisOrientaion(target)
            (ϕ,θ,0.0)
        end
    α,β,γ_ = obtain_FF_MF_Euler(target_rot, [Fxtr,Fytr])
    # computes wigner-D beforehand
    if target isa MoleculeBase
        mol_wigner_D_mat = zeros(ComplexF64, lMax+1, 2lMax+1, 2lMax+1) # to get D_{m,m'}^l (α,β,γ), call wigner_D_mat[l+1, m+l+1, m_+l+1].
        for l in 0:lMax for m in -l:l for m_ in -l:l
            mol_wigner_D_mat[l+1, m+l+1, m_+l+1] = wignerDjmn(l,m,m_, α,β,γ_)
        end; end; end
    elseif target isa SAEAtomBase
        atm_wigner_D_mat = zeros(ComplexF64, 2l+1, 2l+1) # to get D_{m,m'}^l (α,β,γ), call wigner_D_mat[m+l+1, m_+l+1].
        for m in -l:l for m_ in -l:l
            atm_wigner_D_mat[m+l+1, m_+l+1] = wignerDjmn(l,m,m_, α,β,γ_)
        end; end
    end
    # computes spherical harmonics beforehand
    SH_func_mat =   # to get Y_{lm}(x,y,z), call SH_func_mat[l+1,m+l+1]
        if target isa MoleculeBase
            gen_sph_harm_funcs(lMax, sqrt(2Ip))
        elseif target isa SAEAtomBase
            gen_sph_harm_funcs(l, sqrt(2Ip))
        end
    new_x_axis, new_y_axis, new_z_axis = obtain_xyz_FF_LF(Fxtr,Fytr)

    # determines ionization amplitude (which contains phase)
    κ = sqrt(2Ip)
    n = Z/κ # n* = Z/κ
    c = 2^(n/2+1) * κ^(2n+1/2) * gamma(n/2+1) # i^(3(1-n)/2) is omitted, trivial
    e0 = 2.71828182845904523
    c_cc = 2^(3n/2+1) * κ^(5n+1/2) * Ftr^(-n) * (1+2γ/e0)^(-n)
    FxdFy_FydFx = Fxtr*dFyt-dFxt*Fytr
    @inline F2eff(kx,ky) = Ftr^2 - (kx*dFxt+ky*dFyt)  # F2eff=F²-k⟂⋅F'
    @inline r_exit(kx,ky,kz) = Ftr/2*(kx^2+ky^2+kz^2+2Ip)/F2eff(kx,ky)
    @inline ti(kx,ky,kd,kz) = sqrt((κ^2+kd^2+kz^2)/F2eff(kx,ky))
    @inline k_ts(kx,ky,kz,ts) = SVector{3,ComplexF64}(kx-1im*imag(ts)*Fxtr+imag(ts)^2/2*dFxt, ky-1im*imag(ts)*Fytr+imag(ts)^2/2*dFyt, kz)
    if target isa MoleculeBase
        lmm_list = [(l,m,m_) for l in 0:lMax for m in -l:l for m_ in -l:l] # this cache improves the efficiency a lot
    elseif target isa SAEAtomBase
        lmm_list = [(l,m,m_) for m in -l:l for m_ in -l:l]
    end
    @inline function pre_numerator((l,m,m_),kts)
        if target isa MoleculeBase
            C_lm = asymp_coeff[l+1, m+l+1]
            (C_lm == 0) && return 0.0
            new_kx = kts ⋅ new_x_axis
            new_ky = kts ⋅ new_y_axis
            new_kz = kts ⋅ new_z_axis
            C_lm * mol_wigner_D_mat[l+1,m_+l+1,m+l+1] * SH_func_mat[l+1,m_+l+1](new_kx, new_ky, new_kz)
        elseif target isa SAEAtomBase
            new_kx = kts ⋅ new_x_axis
            new_ky = kts ⋅ new_y_axis
            new_kz = kts ⋅ new_z_axis
            C * atm_wigner_D_mat[m_+l+1,m+l+1] * SH_func_mat[l+1,m_+l+1](new_kx, new_ky, new_kz)
        end
    end
    pre(kx,ky,kz,ts) = begin
        kts = k_ts(kx,ky,kz,ts)
        c * mapreduce(((l,m,m_),) -> pre_numerator((l,m,m_),kts), +, lmm_list) / ((kx^2+ky^2+kz^2+2Ip)*F2eff(kx,ky))^((n+1)/4)
    end
    pre_cc(kx,ky,kz,ts) = begin
        kts = k_ts(kx,ky,kz,ts)
        c_cc * mapreduce(((l,m,m_),) -> pre_numerator((l,m,m_),kts), +, lmm_list) / ((kx^2+ky^2+kz^2+2Ip)*F2eff(kx,ky))^((n+1)/4)
    end
    jac(kd,kz) = abs(Ftr + sqrt(kd^2+kz^2)/Ftr^2*FxdFy_FydFx)
    step(range) = (maximum(range)-minimum(range))/length(range) # gets the step length of the range
    dkdt = if ! sp.monte_carlo
        kd_step = length(sp.ss_kd_samples)>1 ? step(sp.ss_kd_samples) : 1.0
        kz_step = (sp.dimension==3 && length(sp.ss_kz_samples)>1) ? step(sp.ss_kz_samples) : 1.0
        step(sp.t_samples) * kd_step * kz_step
    else
        kd_span = sp.mc_kd_max>0 ? 2*sp.mc_kd_max : 1.0
        kz_span = sp.mc_kz_max>0 ? 2*sp.mc_kz_max : 1.0
        step(sp.t_samples) * kd_span*kz_span/sp.mc_kt_num
    end
    amplitude::Function =
        if isempty(prefix)
            amp_exp(kx,ky,kd,kz,ti) = sqrt(dkdt) * ADKAmpExp(sqrt(F2eff(kx,ky)),Ip,kd,kz)
        else
            if :Pre in prefix
                if :Jac in prefix
                    amp_pre_jac(kx,ky,kd,kz,ti) = sqrt(dkdt) * ADKAmpExp(sqrt(F2eff(kx,ky)),Ip,kd,kz) * pre(kx,ky,kz,tr+1im*ti) * sqrt(jac(kd,kz))
                else
                    amp_pre(kx,ky,kd,kz,ti) = sqrt(dkdt) * ADKAmpExp(sqrt(F2eff(kx,ky)),Ip,kd,kz) * pre(kx,ky,kz,tr+1im*ti)
                end
            elseif :PreCC in prefix
                if :Jac in prefix
                    amp_precc_jac(kx,ky,kd,kz,ti) = sqrt(dkdt) * ADKAmpExp(sqrt(F2eff(kx,ky)),Ip,kd,kz) * pre_cc(kx,ky,kz,tr+1im*ti) * sqrt(jac(kd,kz))
                else
                    amp_precc(kx,ky,kd,kz,ti) = sqrt(dkdt) * ADKAmpExp(sqrt(F2eff(kx,ky)),Ip,kd,kz) * pre_cc(kx,ky,kz,tr+1im*ti)
                end
            else # [:Jac]
                amp_jac(kx,ky,kd,kz,ti) = sqrt(dkdt) * ADKAmpExp(sqrt(F2eff(kx,ky)),Ip,kd,kz) * sqrt(jac(kd,kz))
            end
        end

    phase_method = sp.phase_method
    dim = if sp.dimension == 3
        (phase_method == :CTMC) ? 8 : 9 # x,y,z,kx,ky,kz,t0,rate[,phase]
    else
        (phase_method == :CTMC) ? 6 : 7 # x,y,kx,ky,t0,rate[,phase]
    end

    sample_count_thread = zeros(Int,nthreads())
    init_thread = if ! sp.monte_carlo
        zeros(Float64, dim, nthreads(), length(sp.ss_kd_samples)*length(sp.ss_kz_samples)) # initial condition (support for multi-threading)
    else
        zeros(Float64, dim, nthreads(), sp.mc_kt_num)
    end

    k0_cutoff_crit = 1e-4

    if ! sp.monte_carlo
        kd_samples = sp.ss_kd_samples
        kz_samples = sp.ss_kz_samples
        kdNum, kzNum = length(kd_samples), length(kz_samples)
        @threads for ikd in 1:kdNum
            for ikz in 1:kzNum
                kd0, kz0 = kd_samples[ikd], kz_samples[ikz]
                if kd0^2+kz0^2 < k0_cutoff_crit^2
                    continue
                end
                kx0 = kd0*-sin(φ)
                ky0 = kd0* cos(φ)
                r0 = r_exit(kx0,ky0,kz0)
                x0 = r0*cos(φ)
                y0 = r0*sin(φ)
                z0 = 0.0
                if F2eff(kx0,ky0) ≤ Ftr^2/9 # cut off is set to F/3 to avoid NaN probabilities.
                    continue
                end
                amp = amplitude(kx0,ky0,kd0,kz0,ti(kx0,ky0,kd0,kz0))
                rate = abs2(amp)
                if isnan(rate) || rate/dkdt < cutoff_limit
                    continue    # discard the sample
                end
                sample_count_thread[threadid()] += 1
                if sp.dimension == 3
                    init_thread[1:8,threadid(),sample_count_thread[threadid()]] = [x0,y0,z0,kx0,ky0,kz0,tr,rate]
                    if phase_method != :CTMC
                        init_thread[9,threadid(),sample_count_thread[threadid()]] = angle(amp)
                    end
                else # dimension == 2
                    init_thread[1:6,threadid(),sample_count_thread[threadid()]] = [x0,y0,kx0,ky0,tr,rate]
                    if phase_method != :CTMC
                        init_thread[7,threadid(),sample_count_thread[threadid()]] = angle(amp)
                    end
                end
            end
        end
    else
        @threads for i in 1:sp.mc_kt_num
            if sp.dimension == 3
                kd0, kz0 = gen_rand_pt_2dsq(sp.mc_kd_max, sp.mc_kz_max)
                if sp.mc_kd_max == 0
                    kd0 = 0.0
                    kz0 = gen_rand_pt_1d(sp.mc_kz_max)
                elseif sp.mc_kz_max == 0
                    kd0 = gen_rand_pt_1d(sp.mc_kd_max)
                    kz0 = 0.0
                end
            else
                kd0 = gen_rand_pt_1d(sp.mc_kd_max)
                kz0 = 0.0
            end
            if kd0^2+kz0^2 < k0_cutoff_crit^2
                continue
            end
            kx0 = kd0*-sin(φ)
            ky0 = kd0* cos(φ)
            r0 = r_exit(kx0,ky0,kz0)
            x0 = r0*cos(φ)
            y0 = r0*sin(φ)
            z0 = 0.0
            if F2eff(kx0,ky0) ≤ Ftr^2/25
                continue
            end
            amp = amplitude(kx0,ky0,kd0,kz0,ti(kx0,ky0,kd0,kz0))
            rate = abs2(amp)
            if isnan(rate) || rate/dkdt < cutoff_limit
                continue    # discard the sample
            end
            sample_count_thread[threadid()] += 1
            if sp.dimension == 3
                init_thread[1:8,threadid(),sample_count_thread[threadid()]] = [x0,y0,z0,kx0,ky0,kz0,tr,rate]
                if phase_method != :CTMC
                    init_thread[9,threadid(),sample_count_thread[threadid()]] = angle(amp)
                end
            else # dimension == 2
                init_thread[1:6,threadid(),sample_count_thread[threadid()]] = [x0,y0,kx0,ky0,tr,rate]
                if phase_method != :CTMC
                    init_thread[7,threadid(),sample_count_thread[threadid()]] = angle(amp)
                end
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