
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
using OrdinaryDiffEq: ODEProblem, remake, solve, Tsit5, EnsembleProblem, EnsembleSerial, ReturnCode
using ProgressMeter: Progress, ProgressUnknown, BarGlyphs, next!, finish!
using Pkg

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
    batch_completion;
    eff_trajs;
    spec_collect;
    classical_prob;
    classical_prob_uncollected;

    file;
    output_fmt;
    output_compress;
    output_path;
    params;
end

"""
Performs a semiclassical trajectory simulation with given parameters.

# Parameters

## Required parameters:
- `init_cond_method`    : Method used to determine the initial conditions of electrons.
    - Candidates: `:ADK`, `:SPA` (SFA-SPA), `:SPANE` (SFA-SPANE) for targets of type `SAEAtomBase` or `MoleculeBase`; `:WFAT` for `MoleculeBase` targets.
- `laser::Laser`        : A `Lasers.Laser` object which stores parameters of the laser field. See the [Lasers](@ref) module for details.
- `target::Target`      : A `Targets.Target` object which stores parameters of the target. See the [Targets](@ref) module for details.
- `dimension = 2|3`     : Dimensionality of simulation.
    - 2D simulation is carried out in the xy plane.
- `sample_t_intv`       : Time interval for sampling initial electrons.
    - Format: `(start,stop)`
    - Unit: pass numerically in **a.u.** or pass as a `Unitful.Quantity`.
- `sample_t_num`        : Number of time samples.
- `traj_t_final`        : Final time of each trajectory simulation
    - Unit: numerically in **a.u.** or pass as a `Unitful.Quantity`.
- `final_p_max`         : Boundaries of final momentum grid. Grid ranges from `-pxMax` to `+pxMax` in the x direction, and the same for y and z directions.
    - Format: `(pxMax,pyMax[,pzMax])`
- `final_p_num`         : Numbers of final momentum grid points. If a value is `1`, electrons will be collected regardless of the momentum on that dimension.
    - Format: `(pxNum,pyNum[,pzNum])`

## Required parameters for step-sampling methods:
- `ss_kd_max`   : Boundary of k⟂ samples (in a.u.). k⟂ ranges from `-ss_kd_max` to `+ss_kd_max`.
- `ss_kd_num`   : Number of k⟂ samples.
- `ss_kz_max`   : [3D only] Boundary of kz samples (in a.u.). kz ranges from `-ss_kz_max` to `+ss_kz_max`.
- `ss_kz_num`   : [3D only] Number of kz samples.

## Required parameters for Monte-Carlo sampling methods:
- `mc_kt_num`   : Number of kt samples in a single time sample.
- `mc_kd_max`   : Boundary of k⟂ (in a.u.). k⟂ ranges from `-mc_kd_max` to `+mc_kd_max`.
- `mc_kz_max`   : [3D only] Boundary of kz (in a.u.). kz ranges from `-mc_kz_max` to `+mc_kz_max`.

## Optional parameters:
- `traj_phase_method`   : Method used to determine classical trajectories' phase.
    - Candidates: `:CTMC` (default), `:QTMC`, and `:SCTS`.
- `traj_rtol`           : Relative error tolerance for solving classical trajectories (default `1e-6`).
- `output_fmt`          : Output file format.
    - Candidates: `:jld2` (JLD2, default) and `:h5` (HDF5).
- `output_compress`     : Determines whether output files are compressed or not (default `true`).
    - Note: For JLD2 output format, compression requires explicit installation of the `CodecZlib` package.
- `output_path`         : Path to output file.
- `sample_cutoff_limit` : Probability cutoff limit for sampled electrons (default `1e-16`). Electrons with probabilities lower than the limit would be discarded.
- `sample_monte_carlo`  : Determines whether Monte-Carlo sampling is used when generating electron samples (default `false`).

## Optional parameter for atomic SFA-SPA, SFA-SPANE and ADK methods:
- `rate_prefix` : Prefix of the exponential term in the ionization rate (default `:Full`).
    - `:Exp` indicates no prefix; `:Pre` and `:PreCC` indicates inclusion of the prefactor with/without Coulomb correction; `:Jac` indicates inclusion of the Jacobian factor which is related to the sampling method; `:Full` is equivalent to `Set([:PreCC,:Jac])`.
        To combine `:Pre` and `:Jac`, pass `Set([:Pre,:Jac])`.

## Optional parameter for target `MoleculeBase`:
- `mol_orbit_ridx`  : Index of selected orbital relative to the HOMO (e.g., `0` indicates HOMO, and `-1` indicates HOMO-1.)
    - For open-shell molecules, according to α/β spins, should be passed in format `(spin, idx)` where for α orbitals spin=`1` and for β orbitals spin=`2`.

## Optional parameter for interface:
- `show_progress`   : Whether to display progress bar (default `true`).

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
    output_compress     ::Bool      = true,
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
        output_fmt, output_compress, output_path,
        mol_orbit_ridx
        )   # pack up all parameters
    job = TrajectorySimulationJob(;kwargs...)
    batch_num = ElectronSamplers.batch_num(job.sp)
    thread_count = Threads.nthreads()
    prog1 = ProgressUnknown(dt=0.2, desc="Launching electrons and collecting...", color = :cyan, spinner = true, enabled = show_progress)
    prog2 = Progress(batch_num; dt=0.2, color = :cyan, barlen = 25, barglyphs = BarGlyphs('[', '●', ['◔', '◑', '◕'], '○', ']'), showspeed = true, offset=1, enabled = show_progress)
    @threads for i in 1:batch_num
        launch_and_collect!(job, i)
        next!(prog1,spinner=raw"-\|/",desc="Launching electrons and collecting ... [batch #$(sum(job.batch_completion))/$batch_num, $(sum(job.eff_trajs)) electrons collected]"); next!(prog2);
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
        output_fmt, output_compress, output_path,
        ss_kd_max, ss_kd_num, ss_kz_max, ss_kz_num,
        mc_kt_num, mc_kd_max, mc_kz_max,
        traj_phase_method, traj_rtol, sample_cutoff_limit, sample_monte_carlo, rate_prefix,
        mol_orbit_ridx = kwargs
    #* check parameters
    # make conversions
    (traj_t_final isa Quantity) && (traj_t_final = (uconvert(u"fs", traj_t_final) |> auconvert).val)
    ((sample_t_intv[1] isa Quantity) || (sample_t_intv[2] isa Quantity)) && (sample_t_intv = ((uconvert(u"fs", sample_t_intv[1])|>auconvert).val, (uconvert(u"fs", sample_t_intv[2])|>auconvert).val))
    # sample time & traj final time check
    @assert (sample_t_intv[1] < sample_t_intv[2]) "[TrajectorySimulationJob] `sample_t_intv` must be a valid interval."
    @assert (traj_t_final > sample_t_intv[2])  "[TrajectorySimulationJob] `traj_t_final` must be greater than the last time point of `sample_t_intv`."
    Fx_func = LaserFx(laser); Fy_func = LaserFy(laser);
    (Fx_func(sample_t_intv[1])^2 + Fy_func(sample_t_intv[1])^2 > 0.005^2) && @warn("[TrajectorySimulationJob] The laser field is not vanishingly small at the beginning of the sampling time interval. Consider extending the sampling time interval to include the beginning of the laser pulse.")
    (Fx_func(sample_t_intv[2])^2 + Fy_func(sample_t_intv[2])^2 > 0.005^2) && @warn("[TrajectorySimulationJob] The laser field is not vanishingly small at the end of the sampling time interval. Consider extending the sampling time interval to include the end of the laser pulse.")
    (Fx_func(traj_t_final)^2 + Fy_func(traj_t_final)^2 > 0.005^2) && @warn("[TrajectorySimulationJob] The laser field is not vanishingly small at the end of the trajectory simulation. Consider delaying `traj_t_final`.")
    # dimension check
    @assert dimension in (2,3) "[TrajectorySimulationJob] `dimension` must be either 2 or 3."
    @assert length(final_p_max)==length(final_p_num)==dimension "[TrajectorySimulationJob] `length(final_p_max)` and `length(final_p_num)` should match `dimension`."
    # check degenerate orbitals && molecule size
    if target isa GenericMolecule
        energy_levels = MolEnergyLevels(target)
        deg_orb = Targets._lookup_degenerate_orbital(target, mol_orbit_ridx)
        if !isempty(deg_orb)
            @info "[TrajectorySimulationJob] Target `$(TargetName(target))` has degenerate orbitals $(join(deg_orb, ", ", " and "))."
        end
        if length(MolAtoms(target)) > 5 && init_cond_method in [:SPA, :SPANE, :ADK, :SFA, :SFAAE, :MOSPA, :MOSPANE, :MOADK, :MOSFA, :MOSFAAE]
            @warn "[TrajectorySimulationJob] Target `$(TargetName(target))` is a large molecule, the SFA-based methods may not produce accurate results and are much slower. Please consider using the WFAT-CTMC method by setting `init_cond_method = :WFAT` and `phase_method = :CTMC`."
        end
    end
    # check the path in the first
    if isfile(output_path)
        @warn "[TrajectorySimulationJob] File `$output_path` already exists, will save at `$(default_filename())`."
        output_path = default_filename()
    end
    if isfile(output_path * "!")
        @warn "[TrajectorySimulationJob] The temp file `$(output_path)!` exists (may be created by a previous session), will save at `$(default_filename())`."
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
    # check output compression
    if output_fmt == :jld2 && output_compress
        CodecZlib_installed = false
        for (k,v::Pkg.API.PackageInfo) in Pkg.dependencies()
            v.is_direct_dep || continue
            if v.name == "CodecZlib"
                CodecZlib_installed = true
                break
            end
        end
        if !CodecZlib_installed
            @warn "[TrajectorySimulationJob] JLD2 output compression is enabled (`output_compress==true`) but `CodecZlib` is not explicitly installed. Will disable compression."
            output_compress = false
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
    px = final_p_num[1]>1 ? collect(range(-final_p_max[1],final_p_max[1], length=final_p_num[1])) : [0.0]
    py = final_p_num[2]>1 ? collect(range(-final_p_max[2],final_p_max[2], length=final_p_num[2])) : [0.0]
    pz = if dimension == 3
        final_p_num[3]>1 ? collect(range(-final_p_max[3],final_p_max[3], length=final_p_num[3])) : [0.0]
    else
        nothing
    end
    # ionization amplitude (spec_collect for temporary cache)
    spec_collect =
        if traj_phase_method == :CTMC
            [zeros(Float64, final_p_num...) for _ in 1:nthreads]
        else
            [zeros(ComplexF64, final_p_num...) for _ in 1:nthreads]
        end
    # classical prob
    classical_prob = zeros(Float64, nthreads)
    classical_prob_uncollected = zeros(Float64, nthreads)
    # batch_completion & eff_trajs
    batch_completion = zeros(Int, batch_num(sp))
    eff_trajs = zeros(Int, batch_num(sp))
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
    return TrajectorySimulationJob(sp,laser,target,dimension,traj_phase_method,traj_t_final,traj_rtol,final_p_max,final_p_num,px,py,pz,nthreads,batch_completion,eff_trajs,spec_collect,classical_prob,classical_prob_uncollected,file,output_fmt,output_compress,output_path,params)
end

function launch_and_collect!(job::TrajectorySimulationJob, batch_id::Integer)
    if job.dimension == 2
        launch_and_collect_2D!(job, batch_id)
    else
        launch_and_collect_3D!(job, batch_id)
    end
end

function launch_and_collect_2D!(job::TrajectorySimulationJob, batch_id::Integer)
    @assert job.dimension == 2
    (job.batch_completion[batch_id] == 1) && return
    init = gen_electron_batch(job.sp, batch_id)
    isnothing(init) && (job.batch_completion[batch_id] = 1 ; return)
    # introduce shorthands to avoid the `job.` prefix.
    dimension, target, laser, traj_phase_method, traj_t_final, traj_rtol, final_p_num, final_p_max =
        job.dimension, job.target, job.laser, job.traj_phase_method, job.traj_t_final, job.traj_rtol, job.final_p_num, job.final_p_max
    spec_collect, classical_prob, classical_prob_uncollected =
        job.spec_collect, job.classical_prob, job.classical_prob_uncollected

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
    threadid = Threads.threadid()
    spec_collect_current = spec_collect[threadid]

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
    sol = solve(ensemble_prob, Tsit5(), EnsembleSerial(), trajectories=batch_size, adaptive=true, reltol=traj_rtol, save_everystep=false)
    # collect and summarize.
    for i in 1:batch_size
        x0,y0,px0,py0 = sol.u[i].u[ 1 ][1:4]
        x, y, px, py  = sol.u[i].u[end][1:4]
        if sol.u[i].retcode != ReturnCode.Success
            warn_num += 1
            if warn_num < max_warn_num
                @warn "[launch_and_collect!] ODE solver failed (retcode `$(sol.u[i].retcode)`), the electron's initial condition is r0=$([x0,y0]), k0=$([px0,py0]), t0=$(sol.u[i].t[1]), rate=$(init[6,i])."
            elseif warn_num == max_warn_num
                @warn "[launch_and_collect!] ODE solver failed (retcode `$(sol.u[i].retcode)`), the electron's condition is r0=$([x0,y0]), k0=$([px0,py0]), t0=$(sol.u[i].t[1]), rate=$(init[6,i]). Similar warnings in the same batch would be suppressed."
            end
            continue
        end
        phase = (traj_phase_method == :CTMC) ? (0.) : (sol.u[i].u[end][5])
        prob = init[6,i]
        E_inf = (px^2+py^2)/2 + targetP(x,y,0.0)
        r_vec = [x, y, 0.0]
        p_vec = [px,py,0.0]
        L_vec = r_vec × p_vec
        L2    = sum(abs2.(L_vec))
        if E_inf > 0 && traj_phase_method == :SCTS
            phase += Ip * init[5,i] # Ip*tr
            # asymptotic Coulomb phase correction term in SCTS
            sqrtb = (2E_inf)^(-0.5)
            g = sqrt(1+2E_inf*L2)
            phase -= nucl_charge*sqrtb*(log(g)+asinh((x*px+y*py)/(g*sqrtb)))
        end
        if E_inf ≥ 0    # finally ionized.
            classical_prob[threadid] += prob
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
            if checkbounds(Bool, spec_collect_current, pxIdx,pyIdx)
                if traj_phase_method == :CTMC
                    spec_collect_current[pxIdx,pyIdx] += prob # prob
                else
                    spec_collect_current[pxIdx,pyIdx] += sqrt(prob)*exp(1im*phase) # sqrt(prob)*phase_factor
                end
            else
                classical_prob_uncollected[threadid] += prob
            end
        else # finally become rydberg.
            # currently does not support Rydberg collecting.
        end
    end
    job.batch_completion[batch_id] = 1
    job.eff_trajs[batch_id] += batch_size
    return
end

function launch_and_collect_3D!(job::TrajectorySimulationJob, batch_id::Integer)
    @assert job.dimension == 3
    (job.batch_completion[batch_id] == 1) && return
    init = gen_electron_batch(job.sp, batch_id)
    isnothing(init) && (job.batch_completion[batch_id] = 1 ; return)
    # introduce shorthands to avoid the `job.` prefix.
    dimension, target, laser, traj_phase_method, traj_t_final, traj_rtol, final_p_num, final_p_max =
        job.dimension, job.target, job.laser, job.traj_phase_method, job.traj_t_final, job.traj_rtol, job.final_p_num, job.final_p_max
    spec_collect, classical_prob, classical_prob_uncollected =
        job.spec_collect, job.classical_prob, job.classical_prob_uncollected

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
    threadid = Threads.threadid()
    spec_collect_current = spec_collect[threadid]

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
    sol = solve(ensemble_prob, Tsit5(), EnsembleSerial(), trajectories=batch_size, adaptive=true, reltol=traj_rtol, save_everystep=false)
    # collect and summarize.
    for i in 1:batch_size
        x0,y0,z0,px0,py0,pz0 = sol.u[i].u[ 1 ][1:6]
        x, y, z, px, py, pz  = sol.u[i].u[end][1:6]
        if sol.u[i].retcode != ReturnCode.Success
            warn_num += 1
            if warn_num < max_warn_num
                @warn "[launch_and_collect!] ODE solver failed (retcode `$(sol.u[i].retcode)`), the electron's initial condition is r0=$([x0,y0,z0]), k0=$([px0,py0,pz0]), t0=$(sol.u[i].t[1]), rate=$(init[8,i])."
            elseif warn_num == max_warn_num
                @warn "[launch_and_collect!] ODE solver failed (retcode `$(sol.u[i].retcode)`), the electron's condition is r0=$([x0,y0,z0]), k0=$([px0,py0,pz0]), t0=$(sol.u[i].t[1]), rate=$(init[8,i]). Similar warnings in the same batch would be suppressed."
            end
            continue
        end
        phase = (traj_phase_method == :CTMC) ? (0.) : (sol.u[i].u[end][7])
        prob = init[8,i]
        E_inf = (px^2+py^2+pz^2)/2 + targetP(x,y,z)
        r_vec = [x, y, z ]
        p_vec = [px,py,pz]
        L_vec = r_vec × p_vec
        L2    = sum(abs2.(L_vec))
        if E_inf > 0 && traj_phase_method == :SCTS
            phase += Ip * init[7,i] # Ip*tr
            # asymptotic Coulomb phase correction term in SCTS
            sqrtb = (2E_inf)^(-0.5)
            g = sqrt(1+2E_inf*L2)
            phase -= nucl_charge*sqrtb*(log(g)+asinh((x*py+y*py+z*pz)/(g*sqrtb)))
        end
        if E_inf ≥ 0    # finally ionized.
            classical_prob[threadid] += prob
            p_inf = sqrt(2E_inf)
            a_vec = p_vec × L_vec - nucl_charge * r_vec ./ norm(r_vec)
            p_inf_vec = (p_inf/(1+p_inf^2*L2)) .* (p_inf .* (L_vec×a_vec) - a_vec)
            if final_p_num[1]>1
                pxIdx = round(Int, (p_inf_vec[1]+final_p_max[1])/(final_p_max[1]/final_p_num[1]*2))
            else
                pxIdx = 1   # we need to unconditionally collect if final_p_max[i]==1, without this step, only electrons with positive p_x component would be collected, the same for p_y and p_z.
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
            if checkbounds(Bool, spec_collect_current, pxIdx,pyIdx,pzIdx)
                if traj_phase_method == :CTMC
                    spec_collect_current[pxIdx,pyIdx,pzIdx] += prob # prob
                else
                    spec_collect_current[pxIdx,pyIdx,pzIdx] += sqrt(prob)*exp(1im*phase) # sqrt(prob)*phase_factor
                end
            else
                classical_prob_uncollected[threadid] += prob
            end
        else # finally become rydberg.
            # currently does not support Rydberg collecting.
        end
    end
    job.batch_completion[batch_id] = 1
    job.eff_trajs[batch_id] += batch_size
end

function write_output!(job::TrajectorySimulationJob)
    if reduce(*,job.batch_completion) != 1
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
    ion_prob_final = sum(job.spec_collect)
    if job.traj_phase_method == :CTMC
        fmt == :jld2 && write(file, "momentum_spec", ion_prob_final; compress=job.output_compress)
        if fmt == :h5
            if job.output_compress
                file["momentum_spec", shuffle=true, deflate=6] = ion_prob_final
            else
                file["momentum_spec"] = ion_prob_final
            end
        end
    else
        fmt == :jld2 && write(file, "momentum_spec", ion_prob_final .|> abs2; compress=job.output_compress)
        if fmt == :h5
            if job.output_compress
                file["momentum_spec", shuffle=true, deflate=6] = ion_prob_final .|> abs2
            else
                file["momentum_spec"] = ion_prob_final .|> abs2
            end
        end
    end
    file["ion_prob"] = job.classical_prob |> sum
    file["ion_prob_uncollected"] = job.classical_prob_uncollected |> sum
    file["num_effective_traj"] = job.eff_trajs |> sum
    close(file)
    mv(job.output_path * "!", job.output_path)
end

function default_filename()
    Y,M,D = yearmonthday(now())
    h,m,s = hour(now()), minute(now()), second(now())
    return "TrajectorySimulation-$(string(Y,pad=4))$(string(M,pad=2))$(string(D,pad=2)).$(string(h,pad=2))$(string(m,pad=2))$(string(s,pad=2))"
end