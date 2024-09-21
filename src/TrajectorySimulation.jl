
using Dates
using Unitful, UnitfulAtomic
using HDF5: h5open
using JLD2: jldopen, write
using YAML  # YAML: write
using OrderedCollections: OrderedDict
using Parameters: @pack!, @unpack
using Base.Threads
using LinearAlgebra: norm, ×
using StaticArrays: SVector, @SVector
using OrdinaryDiffEq: ODEProblem, remake, solve, Tsit5, EnsembleProblem, EnsembleThreads
using ProgressMeter: Progress, ProgressUnknown, BarGlyphs, next!, finish!

mutable struct TrajectorySimulationJob
    sp::ElectronSampler;
    laser::Laser;
    target::Target;
    dimension;
    traj_phase_method;
    traj_t_final;
    traj_rtol;
    final_p_max;
    final_p_num;
    px;
    py;
    pz;

    nthreads;
    sample_step;
    n_eff_traj;
    ion_prob_sum_temp;
    ion_prob_collect;
    classical_prob;
    ion_prob_final;

    file;
    output_fmt;
    output_path;
    params;
end

"""
Performs a semiclassical simulation with given parameters.

# Parameters

## Required params. for all:
- `init_cond_method = :ADK|:SPA|:SPANE|:WFAT`           : Method of electrons' initial conditions. Currently supports `:ADK`, `:SPA` (SFA-SPA), `:SPANE` (SFA-SPANE) for atoms and molecules, and `:WFAT` for molecules only.
- `laser::Laser`                                        : Parameters of the laser field.
- `target::Target`                                      : Parameters of the target.
- `dimension = 2|3`                                     : Dimension of simulation which indicates 2D/3D simulation, 2D simulation is carried out in the xy plane.
- `sample_t_intv = (start,stop)`                        : Time interval in which the initial electrons are sampled (numerically in **a.u.** or `Unitful.Quantity`).
- `sample_t_num`                                        : Number of time samples.
- `traj_t_final`                                        : Time when every trajectory simulation ends (numerically in **a.u.** or a `Unitful.Quantity`).
- `final_p_max = (pxMax,pyMax[,pzMax])`                 : Boundaries of final momentum spectrum collected in two/three dimensions.
- `final_p_num = (pxNum,pyNum[,pzNum])`                 : Numbers of final momentum spectrum collected in two/three dimensions.

## Required params. for step-sampling methods:
- `ss_kd_max`   : Boundary of kd (momentum's transversal component in the polarization (xy) plane) samples (in a.u.).
- `ss_kd_num`   : Number of kd (momentum's transversal component in the polarization (xy) plane) samples (in a.u.).
- `ss_kz_max`   : [3D simulation] Boundary of kz (momentum's component along propagation direction (z ax.)) samples (in a.u.).
- `ss_kz_num`   : [3D simulation] Number of kz (momentum's component along propagation direction (z ax.)) samples (an even number is required) (in a.u.).

## Required params. for Monte-Carlo-sampling methods:
- `mc_kt_num`   : Number of kt (initial momentum which is perpendicular to field direction, two dimensional) samples in a single time sample.
- `mc_kd_max`   : Boundary of kd.
- `mc_kz_max`   : [3D simulation] Boundary of kz.

## Optional params. for all:
- `traj_phase_method = :CTMC|:QTMC|:SCTS`           : Method of classical trajectories' phase (default `CTMC`).
- `traj_rtol = 1e-6`                                : Relative error tolerance when solving classical trajectories using adaptive methods (default `1e-6`).
- `output_fmt = :jld2|:h5`                          : Output format, two options are JLD2 (`:jld2`) and HDF5 (`:h5`) (Default `:jld2`).
- `output_path`                                     : Output file path.
- `sample_cutoff_limit = 1e-16`                     : The cut-off limit of the probability of the sampled electron, electrons with probabilities lower than the limit would be discarded (Default `1e-16`).
- `sample_monte_carlo = false`                      : Determines whether Monte-Carlo sampling is used when generating electron samples (default `false`).

## Optional params. for atomic SFA-SPA, SFA-SPANE and ADK methods:
- `rate_prefix = :Full | Set([:Pre|:PreCC,:Jac]) | :Exp`    : Prefix of the exponential term in the ionization rate (default `:Full`).
                                                              `:Exp` indicates no prefix; `:Pre` and `:PreCC` indicates inclusion of the prefactor with/without Coulomb correction; `:Jac` indicates inclusion of the Jacobian factor which is related to the sampling method; `:Full` is equivalent to `Set([:PreCC,:Jac])`.
                                                              To combine the prefactor and Jacobian factor, pass a `Set` containing `:Pre` OR `:PreCC`, as well as `:Jac`.

## Optional params. for target `Molecule`:
- `mol_orbit_ridx = 0`  : Index of selected orbit relative to the HOMO (e.g., `0` indicates HOMO, and `-1` indicates HOMO-1).
                          For open-shell molecules, according to α/β spins, should be passed in format `(spin, idx)` where for α orbitals spin=`1` and for β orbitals spin=`2`.

## Optional params. for interface:
- `show_progress = true`    : Indicates whether to show the progressbar (Default `true`).

"""
function perform_traj_simulation(;
    # some abbrs.:  req. = required, opt. = optional, params. = parameters.
        #* req. params. for all methods
    init_cond_method    ::Symbol,
    laser               ::Laser,
    target              ::Target,
    dimension           ::Integer,
    sample_t_intv,      # Tuple{<:Real,<:Real}
    sample_t_num        ::Integer,
    traj_t_final,       # Real
    final_p_max,        # Tuple{<:Real,<:Real} or Tuple{<:Real,<:Real,<:Real}
    final_p_num,        # Tuple{<:Integer,<:Integer} or Tuple{<:Integer,<:Integer,<:Integer}
        #* req. params. for step-sampling (ss) methods
    ss_kd_max           ::Real      = 0.,
    ss_kd_num           ::Integer   = 0 ,
    ss_kz_max           ::Real      = 0.,
    ss_kz_num           ::Integer   = 0 ,
        #* req. params. for Monte-Carlo (mc) methods
    mc_kt_num           ::Integer   = 0 ,
    mc_kd_max           ::Real      = 0.,
    mc_kz_max           ::Real      = 0.,
        #* opt. params. for all methods
    traj_phase_method   ::Symbol    = :CTMC,
    traj_rtol           ::Real      = 1e-6,
    output_fmt          ::Symbol    = :jld2,
    output_path         ::String    = default_filename(),
    sample_monte_carlo  ::Bool      = false,
    sample_cutoff_limit ::Real      = 1e-16,
        #* opt. params. for SFA, SFA-AE and ADK methods
    rate_prefix         ::Union{Symbol,AbstractVector{Symbol},AbstractSet{Symbol}} = :Full,
        #* opt. params. for target `MoleculeBase`
    mol_orbit_ridx                  = 0,
        #* params. for interface
    show_progress       ::Bool      = true
    )
    kwargs = Dict{Symbol,Any}()
    @pack! kwargs= (
        init_cond_method, laser, target, dimension,
        sample_t_intv, sample_t_num, traj_t_final, final_p_max, final_p_num,
        ss_kd_max, ss_kd_num, ss_kz_max, ss_kz_num,
        mc_kt_num, mc_kd_max, mc_kz_max,
        traj_phase_method, traj_rtol, sample_cutoff_limit, sample_monte_carlo, rate_prefix,
        output_fmt, output_path,
        mol_orbit_ridx
        )   # pack up all parameters
    job = TrajectorySimulationJob(;kwargs...)
    batch_num = ElectronSamplers.batch_num(job.sp)
    prog1 = ProgressUnknown(dt=0.2, desc="Launching electrons and collecting...", color = :cyan, spinner = true, enabled = show_progress)
    prog2 = Progress(batch_num; dt=0.2, color = :cyan, barlen = 25, barglyphs = BarGlyphs('[', '●', ['◔', '◑', '◕'], '○', ']'), showspeed = true, offset=1, enabled = show_progress)
    for i in 1:batch_num
        launch_and_collect!(job)
        next!(prog1,spinner=raw"-\|/",desc="Launching electrons and collecting ... [batch #$i/$batch_num, $(job.n_eff_traj) electrons collected]"); next!(prog2);
    end
    finish!(prog1); finish!(prog2); println();
    write_output!(job)
    @info "Task finished, data saved at `$(job.output_path)`."
end

"Initializes a new instance of `TrajectorySimulationJob`."
function TrajectorySimulationJob(; kwargs...)
    #* unpack from kwargs
    @unpack init_cond_method, laser, target, dimension,
        sample_t_intv, sample_t_num, traj_t_final, final_p_max, final_p_num,
        output_fmt, output_path,
        ss_kd_max, ss_kd_num, ss_kz_max, ss_kz_num,
        mc_kt_num, mc_kd_max, mc_kz_max,
        traj_phase_method, traj_rtol, sample_cutoff_limit, sample_monte_carlo, rate_prefix,
        mol_orbit_ridx = kwargs
    #* check parameters
    # make conversions
    (traj_t_final isa Quantity) && (traj_t_final = (uconvert(u"fs", traj_t_final) |> auconvert).val)
    ((sample_t_intv[1] isa Quantity) || (sample_t_intv[2] isa Quantity)) && (sample_t_intv = ((uconvert(u"fs", sample_t_intv[1])|>auconvert).val, (uconvert(u"fs", sample_t_intv[2])|>auconvert).val))
    # dimension check
    @assert dimension in (2,3) "[TrajectorySimulationJob] `dimension` must be either 2 or 3."
    @assert length(final_p_max)==length(final_p_num)==dimension "[TrajectorySimulationJob] `length(final_p_max)` and `length(final_p_num)` should match `dimension`."
    # 3D kz=0 problem
    if dimension == 3 && !sample_monte_carlo && isodd(ss_kz_num)
        @warn "[TrajectorySimulationJob] `ss_kz_num`=$ss_kz_num is an odd number, which may result in anomalous final electron states. Please choose an even number to avoid such problem."
    end
    # check degenerate orbitals
    if target isa GenericMolecule
        energy_levels = MolEnergyLevels(target)
        deg_orb = Targets._lookup_degenerate_orbital(target, mol_orbit_ridx)
        if !isempty(deg_orb)
            @info "[TrajectorySimulationJob] Target `$(TargetName(target))` has degenerate orbitals $(join(deg_orb, ", ", " and "))."
        end
    end

    # check the path in the first
    if isfile(output_path)
        @warn "[TrajectorySimulationJob] File `$output_path` already exists, will save at `$(default_filename())`."
        output_path = default_filename()
    end
    if output_fmt in (:h5, :hdf5)
        output_fmt = :h5
        if !endswith(output_path,".h5")
            output_path *= ".h5"
        end
    elseif output_fmt == :jld2
        if !endswith(output_path,".jld2")
            output_path *= ".jld2"
        end
    end
    # create the file, try to write (and lock it)
    file = if output_fmt == :h5
        h5open(output_path * "!", "w")
    elseif output_fmt == :jld2
        jldopen(output_path * "!", "w")
    else
        error("[TrajectorySimulationJob] Unsupported output format `$output_fmt`.")
    end
    file["info"] = "Trajectory simulation output generated by eTraj @ $(get_version())"
    #* initialize sampler.
    kwargs2 = Dict{Symbol,Any}()
    @pack! kwargs2= (
        init_cond_method, laser, target, dimension,
        sample_t_intv, sample_t_num, traj_t_final, final_p_max, final_p_num,
        ss_kd_max, ss_kd_num, ss_kz_max, ss_kz_num,
        mc_kt_num, mc_kd_max, mc_kz_max,
        traj_phase_method, traj_rtol, sample_cutoff_limit, sample_monte_carlo, rate_prefix,
        mol_orbit_ridx
        )   # pack up all parameters
    sp::ElectronSampler = init_sampler(;kwargs2...)
    #* prepare storage
    nthreads = Threads.nthreads()
    px = collect(range(-final_p_max[1],final_p_max[1], length=final_p_num[1]))
    py = collect(range(-final_p_max[2],final_p_max[2], length=final_p_num[2]))
    pz = if dimension == 3
        collect(range(-final_p_max[2],final_p_max[2], length=final_p_num[2]))
    else
        nothing
    end
    # ionization amplitude (ion_prob_final is for final data, ion_prob_sum_temp & ion_prob_collect are for temporary cache)
    ion_prob_final, ion_prob_sum_temp, ion_prob_collect =
        if traj_phase_method == :CTMC
            zeros(Float64, final_p_num),     zeros(Float64, final_p_num),     zeros(Float64, tuple(final_p_num...,nthreads))
        else
            zeros(ComplexF64, final_p_num),  zeros(ComplexF64, final_p_num),  zeros(ComplexF64, tuple(final_p_num...,nthreads))
        end
    # classical prob
    classical_prob = Dict{Symbol,Float64}()
        classical_prob[:ion]                = 0.
        classical_prob[:ion_uncollected]    = 0.
    #* pack up effective parameters
    params = OrderedDict{Symbol,Any}()
    @pack! params = (
        init_cond_method, laser, target,
        sample_t_intv, sample_t_num, sample_monte_carlo, sample_cutoff_limit,
        traj_phase_method, traj_rtol, traj_t_final,
        output_fmt, output_path,
        final_p_max, final_p_num
    )
    if !sample_monte_carlo
        @pack! params = (ss_kd_max,ss_kd_num)
        dimension==3 && @pack! params = (ss_kz_max,ss_kz_num)
    else
        @pack! params = (mc_kd_max,mc_kt_num)
        dimension==3 && @pack! params = (mc_kz_max)
    end
    if mapreduce(str->endswith(string(init_cond_method), str), |, ("ADK","SPA","SPANE","SFA","SFAAE"))
        @pack! params = (rate_prefix)
    end
    if target isa MoleculeBase
        @pack! params = (mol_orbit_ridx)
    end
    #* initialize job
    return TrajectorySimulationJob(sp,laser,target,dimension,traj_phase_method,traj_t_final,traj_rtol,final_p_max,final_p_num,px,py,pz,nthreads,0,0,ion_prob_sum_temp,ion_prob_collect,classical_prob,ion_prob_final,file,output_fmt,output_path,params)
end

function launch_and_collect!(job::TrajectorySimulationJob)
    if job.dimension == 2
        launch_and_collect_2D!(job)
    else
        launch_and_collect_3D!(job)
    end
end

function launch_and_collect_2D!(job::TrajectorySimulationJob)
    @assert job.dimension == 2
    if job.sample_step >= batch_num(job.sp)
        @warn "[launch_and_collect!] `TrajectorySimulationJob` has already finished collection of all $(job.sample_step) batches of electrons, no electrons would be launched."
        return
    end
    init = gen_electron_batch(job.sp, job.sample_step+1)
    isnothing(init) && (job.sample_step += 1; return)
    # introduce shorthands to avoid the `job.` prefix.
    dimension, target, laser, traj_phase_method, traj_t_final, traj_rtol, final_p_num, final_p_max =
        job.dimension, job.target, job.laser, job.traj_phase_method, job.traj_t_final, job.traj_rtol, job.final_p_num, job.final_p_max
    ion_prob_sum_temp, ion_prob_collect, classical_prob, ion_prob_final =
        job.ion_prob_sum_temp, job.ion_prob_collect, job.classical_prob, job.ion_prob_final

    Ip                  = IonPotential(target)
    Fx::Function        = LaserFx(laser)
    Fy::Function        = LaserFy(laser)
    targetF::Function   = TargetForce(target)
    targetP::Function   = TargetPotential(target)
    nucl_charge         = AsympNuclCharge(target)
    traj::Function      = TrajectoryFunction(target, dimension, Fx,Fy,traj_phase_method)
    batch_size = size(init,2)
    warn_num = 0    # number of warnings of anomalous electrons.
    max_warn_num = 5
    nthreads = Threads.nthreads()
    class_prob_ion              = zeros(nthreads)
    class_prob_ion_uncollected  = zeros(nthreads)

    # create ODE problem and solve the ensemble.
    prob_dim = (traj_phase_method == :CTMC) ? 4 : 5 # x,y,px,py[,phase]
    traj_ODE_prob::ODEProblem = ODEProblem(traj, (@SVector zeros(Float64,prob_dim)), Float64.((0,traj_t_final)))
    init_traj::Function =
        if prob_dim == 4
            (prob,i,repeat) -> remake(prob; u0=SVector{4}([init[k,i] for k in 1:4]),     tspan = (init[5,i],Float64(traj_t_final)))
        else
            (prob,i,repeat) -> remake(prob; u0=SVector{5}([init[k,i] for k in [1:4;7]]), tspan = (init[5,i],Float64(traj_t_final)))
        end
    ensemble_prob = EnsembleProblem(traj_ODE_prob, prob_func=init_traj, safetycopy=false)
    sol = solve(ensemble_prob, Tsit5(), EnsembleThreads(), trajectories=batch_size, adaptive=true, dt=0.01, reltol=traj_rtol, save_everystep=false)
    # collect and summarize.
    Threads.@threads for i in 1:batch_size
        threadid = Threads.threadid()
        x0,y0,px0,py0 = sol.u[i].u[ 1 ][1:4]
        x, y, px, py  = sol.u[i].u[end][1:4]
        if px^2+py^2>100  # possibly anomalous electron, intercept and cancel.
            warn_num += 1
            if warn_num < max_warn_num
                @warn "[launch_and_collect!] Found electron with anomalously large momentum $([px,py]), whose initial condition is r0=$([x0,y0]), k0=$([px0,py0]), t0=$(sol.u[i].t[1])."
            elseif warn_num == max_warn_num
                @warn "[launch_and_collect!] Found electron with anomalously large momentum $([px,py]), whose initial condition is r0=$([x0,y0]), k0=$([px0,py0]), t0=$(sol.u[i].t[1]). Similar warnings in the same batch would be suppressed."
            end
            continue
        end
        phase = (traj_phase_method == :CTMC) ? (0.) : (sol.u[i].u[end][5])
        prob = init[6,i]
        if traj_phase_method == :SCTS # asymptotic Coulomb phase correction term in SCTS
            sqrtb = (2Ip)^(-0.5)
            g = sqrt(1+2Ip*(x*py-y*px)^2)
            phase -= px0*x0+py0*y0 + nucl_charge*sqrtb*(log(g)+asinh((x*px+y*py)/(g*sqrtb)))
        end
        E_inf = (px^2+py^2)/2 + targetP(x,y,0.0)
        r_vec = [x, y, 0.0]
        p_vec = [px,py,0.0]
        L_vec = r_vec × p_vec
        L2    = sum(abs2.(L_vec))
        if E_inf ≥ 0    # finally ionized.
            class_prob_ion[threadid] += prob
            p_inf = sqrt(2E_inf)
            a_vec = p_vec × L_vec - nucl_charge * r_vec ./ norm(r_vec)
            p_inf_vec = (p_inf/(1+p_inf^2*L2)) .* (p_inf .* (L_vec×a_vec) - a_vec)
            if final_p_num[1]>1
                pxIdx = round(Int, (p_inf_vec[1]+final_p_max[1])/(final_p_max[1]/final_p_num[1]*2))
            else
                pxIdx = 1   # without this step, only electrons with positive p_x component would be collected, the same for p_y and p_z.
            end
            if final_p_num[2]>1
                pyIdx = round(Int, (p_inf_vec[2]+final_p_max[2])/(final_p_max[2]/final_p_num[2]*2))
            else
                pyIdx = 1
            end
            if checkbounds(Bool, ion_prob_collect, pxIdx,pyIdx, threadid)
                if traj_phase_method == :CTMC
                    ion_prob_collect[pxIdx,pyIdx, threadid] += prob # prob
                else
                    ion_prob_collect[pxIdx,pyIdx, threadid] += sqrt(prob)*exp(1im*phase) # sqrt(prob)*phase_factor
                end
            else
                class_prob_ion_uncollected[threadid] += prob
            end
        else # finally become rydberg.
            # currently does not support Rydberg collecting.
        end
    end
    sum!(ion_prob_sum_temp, ion_prob_collect)
    ion_prob_final .+= ion_prob_sum_temp
    classical_prob[:ion]                += sum(class_prob_ion)
    classical_prob[:ion_uncollected]    += sum(class_prob_ion_uncollected)

    job.sample_step += 1
    job.n_eff_traj += size(init, 2)
end

function launch_and_collect_3D!(job::TrajectorySimulationJob)
    @assert job.dimension == 3
    if job.sample_step >= batch_num(job.sp)
        @warn "[launch_and_collect!] `TrajectorySimulationJob` has already finished collection of all $(job.sample_step) batches of electrons, no electrons would be launched."
        return
    end
    init = gen_electron_batch(job.sp, job.sample_step+1)
    isnothing(init) && (job.sample_step += 1; return)
    # introduce shorthands to avoid the `job.` prefix.
    dimension, target, laser, traj_phase_method, traj_t_final, traj_rtol, final_p_num, final_p_max =
        job.dimension, job.target, job.laser, job.traj_phase_method, job.traj_t_final, job.traj_rtol, job.final_p_num, job.final_p_max
    ion_prob_sum_temp, ion_prob_collect, classical_prob, ion_prob_final =
        job.ion_prob_sum_temp, job.ion_prob_collect, job.classical_prob, job.ion_prob_final

    Ip                  = IonPotential(target)
    Fx::Function        = LaserFx(laser)
    Fy::Function        = LaserFy(laser)
    targetF::Function   = TargetForce(target)
    targetP::Function   = TargetPotential(target)
    nucl_charge         = AsympNuclCharge(target)
    traj::Function      = TrajectoryFunction(target, dimension, Fx,Fy,traj_phase_method)
    batch_size = size(init,2)
    warn_num = 0    # number of warnings of anomalous electrons.
    max_warn_num = 5
    nthreads = Threads.nthreads()
    class_prob_ion              = zeros(nthreads)
    class_prob_ion_uncollected  = zeros(nthreads)

    # create ODE problem and solve the ensemble.
    prob_dim = (traj_phase_method == :CTMC) ? 6 : 7 # x,y,z,px,py,pz[,phase]
    traj_ODE_prob::ODEProblem = ODEProblem(traj, (@SVector zeros(Float64,prob_dim)), Float64.((0,traj_t_final)))
    init_traj::Function =
        if prob_dim == 6
            (prob,i,repeat) -> remake(prob; u0=SVector{6}([init[k,i] for k in 1:6]),     tspan = (init[7,i],Float64(traj_t_final)))
        else
            (prob,i,repeat) -> remake(prob; u0=SVector{7}([init[k,i] for k in [1:6;9]]), tspan = (init[7,i],Float64(traj_t_final)))
        end
    ensemble_prob = EnsembleProblem(traj_ODE_prob, prob_func=init_traj, safetycopy=false)
    sol = solve(ensemble_prob, Tsit5(), EnsembleThreads(), trajectories=batch_size, adaptive=true, dt=0.01, reltol=traj_rtol, save_everystep=false)
    # collect and summarize.
    Threads.@threads for i in 1:batch_size
        threadid = Threads.threadid()
        x0,y0,z0,px0,py0,pz0 = sol.u[i].u[ 1 ][1:6]
        x, y, z, px, py, pz  = sol.u[i].u[end][1:6]
        if px^2+py^2+pz^2>100  # possibly anomalous electron, intercept and cancel.
            warn_num += 1
            if warn_num < max_warn_num
                @warn "[launch_and_collect!] Found electron with anomalously large momentum $([px,py,pz]), whose initial condition is r0=$([x0,y0,z0]), k0=$([px0,py0,pz0]), t0=$(sol.u[i].t[1])."
            elseif warn_num == max_warn_num
                @warn "[launch_and_collect!] Found electron with anomalously large momentum $([px,py,pz]), whose initial condition is r0=$([x0,y0,z0]), k0=$([px0,py0,pz0]), t0=$(sol.u[i].t[1]). Similar warnings in the same batch would be suppressed."
            end
            continue
        end
        phase = (traj_phase_method == :CTMC) ? (0.) : (sol.u[i].u[end][7])
        prob = init[8,i]
        if traj_phase_method == :SCTS # asymptotic Coulomb phase correction term in SCTS
            sqrtb = (2Ip)^(-0.5)
            g = sqrt(1+2Ip*((y*pz-z*py)^2+(z*px-x*pz)^2+(x*py-y*px)^2))
            phase -= px0*x0+py0*y0+pz0*z0 + nucl_charge*sqrtb*(log(g)+asinh((x*py+y*py+z*pz)/(g*sqrtb)))
        end
        E_inf = (px^2+py^2+pz^2)/2 + targetP(x,y,z)
        r_vec = [x, y, z ]
        p_vec = [px,py,pz]
        L_vec = r_vec × p_vec
        L2    = sum(abs2.(L_vec))
        if E_inf ≥ 0    # finally ionized.
            class_prob_ion[threadid] += prob
            p_inf = sqrt(2E_inf)
            a_vec = p_vec × L_vec - nucl_charge * r_vec ./ norm(r_vec)
            p_inf_vec = (p_inf/(1+p_inf^2*L2)) .* (p_inf .* (L_vec×a_vec) - a_vec)
            if final_p_num[1]>1
                pxIdx = round(Int, (p_inf_vec[1]+final_p_max[1])/(final_p_max[1]/final_p_num[1]*2))
            else
                pxIdx = 1   # without this step, only electrons with positive p_x component would be collected, the same for p_y and p_z.
            end
            if final_p_num[2]>1
                pyIdx = round(Int, (p_inf_vec[2]+final_p_max[2])/(final_p_max[2]/final_p_num[2]*2))
            else
                pyIdx = 1
            end
            if final_p_num[3]>1
                pzIdx = round(Int, (p_inf_vec[3]+final_p_max[3])/(final_p_max[3]/final_p_num[3]*2))
            else
                pzIdx = 1
            end
            if checkbounds(Bool, ion_prob_collect, pxIdx,pyIdx,pzIdx, threadid)
                if traj_phase_method == :CTMC
                    ion_prob_collect[pxIdx,pyIdx,pzIdx, threadid] += prob # prob
                else
                    ion_prob_collect[pxIdx,pyIdx,pzIdx, threadid] += sqrt(prob)*exp(1im*phase) # sqrt(prob)*phase_factor
                end
            else
                class_prob_ion_uncollected[threadid] += prob
            end
        else # finally become rydberg.
            # currently does not support Rydberg collecting.
        end
    end
    sum!(ion_prob_sum_temp, ion_prob_collect)
    ion_prob_final .+= ion_prob_sum_temp
    classical_prob[:ion]                += sum(class_prob_ion)
    classical_prob[:ion_uncollected]    += sum(class_prob_ion_uncollected)

    job.sample_step += 1
    job.n_eff_traj += size(init, 2)
end

function write_output!(job::TrajectorySimulationJob)
    if job.sample_step != batch_num(job.sp)
        @error "[write_output!] The job is not finished yet."
        return
    end
    fmt = job.output_fmt
    file = job.file
    params = copy(job.params)   # the laser, target are serialized into `Dict`s.
    params[:laser] = Lasers.Serialize(params[:laser])
    params[:target] = Targets.Serialize(params[:target])
    file["params_text"] = YAML.write(params) # YAML.write
    if fmt == :jld2
        file["params"] = job.params # h5 does not support saving julia data format
    end
    file["px"] = job.px
    file["py"] = job.py
    if job.dimension == 3
        file["pz"] = job.pz
    end
    if job.traj_phase_method == :CTMC
        fmt == :jld2 && write(file, "momentum_spec", job.ion_prob_final; compress=true)
        fmt == :h5 && (file["momentum_spec", shuffle=true, deflate=6] = job.ion_prob_final)
    else
        fmt == :jld2 && write(file, "momentum_spec", job.ion_prob_final .|> abs2; compress=true)
        fmt == :h5 && (file["momentum_spec", shuffle=true, deflate=6] = job.ion_prob_final .|> abs2)
    end
    file["ion_prob"] = job.classical_prob[:ion]
    file["ion_prob_uncollected"] = job.classical_prob[:ion_uncollected]
    file["num_effective_traj"] = job.n_eff_traj
    close(file)
    mv(job.output_path * "!", job.output_path)
end

function default_filename()
    Y,M,D = yearmonthday(now())
    h,m,s = hour(now()), minute(now()), second(now())
    return "TrajectorySimulation-$(string(Y,pad=4))$(string(M,pad=2))$(string(D,pad=2)).$(string(h,pad=2))$(string(m,pad=2))$(string(s,pad=2))"
end