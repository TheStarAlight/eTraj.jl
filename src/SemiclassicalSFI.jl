
"""
Implementation of semiclassical methods in strong field ionization.
"""
module SemiclassicalSFI

using OrdinaryDiffEq
using DiffEqGPU, CUDA
using LinearAlgebra
using StaticArrays
using Parameters
using HDF5
using Dates
using ProgressMeter
using YAML, OrderedCollections
using Pkg

include("Lasers/Lasers.jl")
include("Targets/Targets.jl")
include("SampleProviders/SampleProviders.jl")
using .Lasers
using .Targets
using .SampleProviders

export performSFI, Lasers, Targets

"""
Performs a semiclassical simulation with given parameters.

# Parameters

## Required params. for all methods:
- `init_cond_method = <:ADK|:SFA|:SFAAE|:WFAT|:MOADK>`  : Method of electrons' initial conditions. Currently supports `:ADK`, `:SFA`, `:SFAAE` for atoms and `:WFAT`, `:MOADK` for molecules.
- `laser::Laser`                                        : Parameters of the laser field.
- `target::Target`                                      : Parameters of the target.
- `sample_t_interval = (start,stop)`                    : Time interval in which the initial electrons are sampled.
- `sample_t_num`                                        : Number of time samples.
- `traj_t_final`                                        : Time when every trajectory simulation ends.
- `final_p_max = (pxMax,pyMax,pzMax)`                   : Boundaries of final momentum spectrum collected in three dimensions.
- `final_p_num = (pxNum,pyNum,pzNum)`                   : Numbers of final momentum spectrum collected in three dimensions.

## Required params. for step-sampling methods:
- `ss_kd_max`   : Boundary of kd (momentum's component along transverse direction (in xy plane)) samples.
- `ss_kd_num`   : Number of kd (momentum's component along transverse direction (in xy plane)) samples.
- `ss_kz_max`   : Boundary of kz (momentum's component along propagation direction (z ax.)) samples.
- `ss_kz_num`   : Number of kz (momentum's component along propagation direction (z ax.)) samples.

## Required params. for Monte-Carlo-sampling methods:
- `mc_kp_num`   : Number of kp (initial momentum which is perpendicular to field direction, two dimensional) samples in a single time sample.
- `mc_kp_max`   : Maximum value of momentum's transversal component (perpendicular to field direction).

## Optional params. for all methods:
- `save_path`                                       : Output HDF5 file path.
- `save_3D_spec = false`                            : Determines whether the 3D momentum spectrum is saved (if not, will save 2D) (default `false`).
- `traj_phase_method = <:CTMC|:QTMC|:SCTS>`         : Method of classical trajectories' phase (default `CTMC`). Currently `:QTMC` and `:SCTS` only supports atom targets.
- `traj_dt = 0.1`                                   : Time step when solving classical trajectories (default `0.1`).
- `traj_nondipole = false`                          : Determines whether non-dipole effect is taken account in the simulation (default `false`).
- `traj_GPU = false`                                : [Experimental] Determines whether GPU acceleration in trajectory simulation is used (default `false`).
- `sample_monte_carlo = false`                      : Determines whether Monte-Carlo sampling is used when generating electron samples (default `false`). Currently only supports ADK.
- `final_ryd_collect = false`                       : Determines whether rydberg final states are collected (default `false`).
- `final_ryd_n_max`                                 : Determines the maximum principle quantum number n for rydberg final states to be collected.

## Optional params. for atomic SFA, SFA-AE and ADK methods:
- `rate_prefix = <:ExpRate|:ExpPre|:ExpJac|:Full>`  : Prefix of the exponential term in the ionization rate (default `:ExpRate`).

## Optional params. for target `Molecule`:
- `mol_orbit_idx = 0`   : Index of the ionizing orbit relative to the HOMO (e.g., `0` indicates HOMO, and `-1` indicates HOMO-1) (default `0`).

## Optional params. for MO-ADK method:
- `moadk_orbit_m = 0`   : Magnetic quantum number m of the ionizing orbital along the z axis. m = 0,1,2 indicate σ, π and δ respectively (default `0`).

## Optional params. for ADK method:
- `adk_tun_exit = <:IpF|:FDM|:Para>` : Tunneling exit method for ADK methods (when `init_cond_method==:ADK`) (default `:IpF`).

"""
function performSFI(; # some abbrs.:  req. = required, opt. = optional, params. = parameters.
                        #* req. params. for all methods
                    init_cond_method    ::Symbol,
                    laser               ::Laser,
                    target              ::Target,
                    sample_t_intv       ::Tuple{<:Real,<:Real},
                    sample_t_num        ::Integer,
                    traj_t_final        ::Real,
                    final_p_max         ::Tuple{<:Real,<:Real,<:Real},
                    final_p_num         ::Tuple{<:Int,<:Int,<:Int},
                        #* req. params. for step-sampling (ss) methods
                    ss_kd_max           ::Real      = 0.,
                    ss_kd_num           ::Integer   = 0 ,
                    ss_kz_max           ::Real      = 0.,
                    ss_kz_num           ::Integer   = 0 ,
                        #* req. params. for Monte-Carlo (mc) methods
                    mc_kp_num           ::Integer   = 0 ,
                    mc_kp_max           ::Real      = 0.,
                        #* opt. params. for all methods
                    save_path           ::String    = default_filename(),
                    save_3D_spec        ::Bool      = false,
                    traj_phase_method   ::Symbol    = :CTMC,
                    traj_dt             ::Real      = 0.1,
                    traj_nondipole      ::Bool      = false,
                    traj_GPU            ::Bool      = false,
                    sample_monte_carlo  ::Bool      = false,
                    final_ryd_collect   ::Bool      = false,
                    final_ryd_n_max     ::Integer   = 0,
                        #* opt. params. for atomic SFA, SFA-AE and ADK methods
                    rate_prefix         ::Symbol    = :ExpRate,
                        #* opt. params. for target `Molecule`
                    mol_orbit_idx       ::Integer   = 0,
                        #* opt. params. for molecular MOADK method
                    moadk_orbit_m       ::Integer   = 0,
                        #* opt. params. for atomic ADK method
                    adk_tun_exit        ::Symbol    = :IpF
                    )
    #* pack up all parameters.
    kwargs = Dict{Symbol,Any}()
    @pack! kwargs= (init_cond_method, laser, target, sample_t_intv, sample_t_num, traj_t_final, final_p_max, final_p_num,
                    ss_kd_max, ss_kd_num, ss_kz_max, ss_kz_num,
                    mc_kp_num, mc_kp_max,
                    traj_phase_method, traj_dt, traj_nondipole, traj_GPU, sample_monte_carlo, rate_prefix, final_ryd_collect, final_ryd_n_max,
                    mol_orbit_idx,
                    moadk_orbit_m,
                    adk_tun_exit)
    #* initialize sample provider.
    sp::ElectronSampleProvider = init_sampler(;kwargs...)
    #* launch electrons and summarize.
    #   * prepare storage
    nthreads = Threads.nthreads()
    # ionization amplitude (ion_prob_final is for final data, ion_prob_sum_temp & ion_prob_collect is for temporary cache)
    ion_prob_final, ion_prob_sum_temp, ion_prob_collect =
        if traj_phase_method == :CTMC
            zeros(Float64, final_p_num),     zeros(Float64, final_p_num),     zeros(Float64, tuple(final_p_num...,nthreads))
        else
            zeros(ComplexF64, final_p_num),  zeros(ComplexF64, final_p_num),  zeros(ComplexF64, tuple(final_p_num...,nthreads))
        end
    # rydberg amplitude
    ryd_prob_final, ryd_prob_sum_temp, ryd_prob_collect =
        if final_ryd_collect
            if traj_phase_method == :CTMC
                zeros(Float64, final_ryd_n_max, final_ryd_n_max, 2*final_ryd_n_max+1),
                zeros(Float64, final_ryd_n_max, final_ryd_n_max, 2*final_ryd_n_max+1),
                zeros(Float64, final_ryd_n_max, final_ryd_n_max, 2*final_ryd_n_max+1, nthreads)
            else
                zeros(ComplexF64, final_ryd_n_max, final_ryd_n_max, 2*final_ryd_n_max+1),
                zeros(ComplexF64, final_ryd_n_max, final_ryd_n_max, 2*final_ryd_n_max+1),
                zeros(ComplexF64, final_ryd_n_max, final_ryd_n_max, 2*final_ryd_n_max+1, nthreads)
            end
        else
            nothing, nothing, nothing
        end
    # classical rates
    classical_rates = Dict{Symbol,Float64}()
        classical_rates[:ion]                = 0.
        classical_rates[:ion_uncollected]    = 0.
        classical_rates[:ryd]                = 0.
        classical_rates[:ryd_uncollected]    = 0.
    #   * launch electrons and collect
    prog1 = ProgressUnknown(dt=0.2, desc="Launching electrons and collecting...", color = :cyan, spinner = true)
    prog2 = Progress(batch_num(sp); dt=0.2, color = :cyan, barlen = 25, barglyphs = BarGlyphs('[', '●', ['◔', '◑', '◕'], '○', ']'), showspeed = true, offset=1)
    warn_num = 0    # number of warnings of empty batches.
    max_warn_num = 5
    for batchId in 1:batch_num(sp)
        init = gen_electron_batch(sp, batchId)
        if isnothing(init)
            warn_num += 1
            if warn_num < max_warn_num
                @warn "[performSFI] The electron sample provider yields no electron sample in batch #$batchId, probably due to zero field strength."
            elseif warn_num == max_warn_num
                @warn "[performSFI] The electron sample provider yields no electron sample in batch #$batchId, probably due to zero field strength. Similar warnings would be suppressed."
            end
        else
            launch_and_collect!(init,
                                ion_prob_final, ion_prob_sum_temp, ion_prob_collect,
                                ryd_prob_final, ryd_prob_sum_temp, ryd_prob_collect,
                                classical_rates; kwargs...)
        end
        next!(prog1,spinner=raw"-\|/"); next!(prog2);
    end
    finish!(prog1); finish!(prog2); println();
    if traj_phase_method != :CTMC
        ion_prob_final = abs2.(ion_prob_final)
        if final_ryd_collect
            ryd_prob_final = abs2.(ryd_prob_final)
        end
    end
    #* save as HDF5.
    if isfile(save_path)
        @warn "[performSFI] File \"$save_path\" already exists. Saving at \"$(default_filename())\"."
        save_path = "$(default_filename()).h5"
    end
    begin
        dict_out = OrderedDict{Symbol,Any}()
        # package version
        dep = Pkg.dependencies()
        for (k,v::Pkg.API.PackageInfo) in dep
            if v.name == "SemiclassicalSFI"
                dict_out[:version] = v.version
            end
        end
        # req. params. for all methods
        dict_out[:init_cond_method]     = init_cond_method
        dict_out[:laser]                = Lasers.Serialize(laser)
        dict_out[:target]               = Targets.Serialize(target)
        dict_out[:sample_t_intv]        = sample_t_intv
        dict_out[:sample_t_num]         = sample_t_num
        dict_out[:traj_t_final]         = traj_t_final
        dict_out[:final_p_max]          = final_p_max
        dict_out[:final_p_num]          = final_p_num
        if ! sample_monte_carlo
            # req. params. for step-sampling (ss) methods
            dict_out[:ss_kd_max]        = ss_kd_max
            dict_out[:ss_kd_num]        = ss_kd_num
            dict_out[:ss_kz_max]        = ss_kz_max
            dict_out[:ss_kz_num]        = ss_kz_num
        else
            # req. params. for Monte-Carlo (mc) methods
            dict_out[:mc_kp_num]        = mc_kp_num
            dict_out[:mc_kp_max]        = mc_kp_max
        end
        # opt. params. for all methods
        dict_out[:save_path]            = save_path
        dict_out[:save_3D_spec]         = save_3D_spec
        dict_out[:traj_phase_method]    = traj_phase_method
        dict_out[:traj_dt]              = traj_dt
        dict_out[:traj_nondipole]       = traj_nondipole
        dict_out[:traj_GPU]             = traj_GPU
        dict_out[:sample_monte_carlo]   = sample_monte_carlo
        dict_out[:final_ryd_collect]    = final_ryd_collect
        dict_out[:final_ryd_n_max]      = final_ryd_n_max
        # opt. params. for atomic SFA, SFA-AE and ADK methods
        if init_cond_method in [:ADK, :SFAAE, :SFA]
            dict_out[:rate_prefix]      = rate_prefix
        end
        # opt. params. for target `Molecule`
        if typeof(target) <: Targets.Molecule
            dict_out[:mol_orbit_idx]    = mol_orbit_idx
        end
        # opt. params. for molecular MOADK method
        if init_cond_method == :MOADK
            dict_out[:moadk_orbit_m]    = moadk_orbit_m
        end
        # opt. params. for atomic ADK method
        if init_cond_method == :ADK
            dict_out[:adk_tun_exit]     = adk_tun_exit
        end
        yaml_out = YAML.write(dict_out)
        h5write(save_path, "abstract", yaml_out)
    end
    h5write(save_path, "px", collect(range(-final_p_max[1],final_p_max[1], length=final_p_num[1])))
    h5write(save_path, "py", collect(range(-final_p_max[2],final_p_max[2], length=final_p_num[2])))
    if save_3D_spec
        h5write(save_path, "pz", collect(range(-final_p_max[3],final_p_max[3], length=final_p_num[3])))
        h5write(save_path, "momentum_spec_3D", ion_prob_final)
    end
    h5write(save_path, "momentum_spec_2D", reshape(sum(ion_prob_final, dims=3),size(ion_prob_final)[1:2]))
    h5write(save_path, "ion_rate",              classical_rates[:ion])
    h5write(save_path, "ion_rate_uncollected",  classical_rates[:ion_uncollected])
    if final_ryd_collect
        h5write(save_path, "ryd_spec", ryd_prob_final)
        h5write(save_path, "ryd_rate",              classical_rates[:ryd])
        h5write(save_path, "ryd_rate_uncollected",  classical_rates[:ryd_uncollected])
    end
    @info "Task finished, data saved at \"$(save_path)\"."
end


function launch_and_collect!( init,
                            ion_prob_final,
                            ion_prob_sum_temp,
                            ion_prob_collect,
                            ryd_prob_final,
                            ryd_prob_sum_temp,
                            ryd_prob_collect,
                            classical_rates;
                            laser               ::Laser,
                            target              ::Target,
                            traj_t_final        ::Real,
                            traj_phase_method   ::Symbol,
                            traj_dt             ::Real,
                            traj_nondipole      ::Bool,
                            traj_GPU            ::Bool,
                            final_ryd_collect   ::Bool,
                            final_ryd_n_max     ::Int,
                            final_p_max         ::Tuple{<:Real,<:Real,<:Real},
                            final_p_num         ::Tuple{<:Int,<:Int,<:Int},
                            kwargs...   # kwargs are surplus params.
                            )
    Ip                  = IonPotential(target)
    Fx::Function        = LaserFx(laser)
    Fy::Function        = LaserFy(laser)
    targetF::Function   = TargetForce(target)
    targetP::Function   = TargetPotential(target)
    nucl_charge         = AsympNuclCharge(target)
    traj::Function      = TrajectoryFunction(target, Fx,Fy,traj_phase_method,traj_nondipole)
    batch_size = size(init,2)
    warn_num = 0    # number of warnings of anomalous electrons.
    max_warn_num = 5
    nthreads = Threads.nthreads()
    class_rates_ion              = zeros(nthreads)
    class_rates_ion_uncollected  = zeros(nthreads)
    class_rates_ryd              = zeros(nthreads)
    class_rates_ryd_uncollected  = zeros(nthreads)

    # create ODE problem and solve the ensemble.
    prob_dim = (traj_phase_method == :CTMC) ? 6 : 7 # x,y,z,px,py,pz[,phase]
    trajODEProb::ODEProblem = ODEProblem(traj, (@SVector zeros(Float64,prob_dim)), (0,traj_t_final))
    initTraj::Function =
        if prob_dim == 6
            (prob,i,repeat) -> remake(prob; u0=SVector{6}([init[k,i] for k in 1:6]),     tspan = (init[7,i],traj_t_final))
        else
            (prob,i,repeat) -> remake(prob; u0=SVector{7}([init[k,i] for k in [1:6;9]]), tspan = (init[7,i],traj_t_final))
        end
    ensemble_prob::EnsembleProblem = EnsembleProblem(trajODEProb, prob_func=initTraj, safetycopy=false)
    sol =
        if ! traj_GPU
            solve(ensemble_prob, OrdinaryDiffEq.Tsit5(), EnsembleThreads(), trajectories=batch_size, adaptive=false, dt=traj_dt, save_everystep=false)
        else
            solve(ensemble_prob, DiffEqGPU.GPUTsit5(), EnsembleGPUKernel(CUDA.CUDABackend()), trajectories=batch_size, adaptive=false, dt=traj_dt, save_everystep=false)
        end
    # collect and summarize.
    Threads.@threads for i in 1:batch_size
        threadid = Threads.threadid()
        x0,y0,z0,px0,py0,pz0 = sol.u[i][ 1 ][1:6]
        x, y, z, px, py, pz  = sol.u[i][end][1:6]
        if px^2+py^2+pz^2>(final_p_max[1]^2+final_p_max[2]^2+final_p_max[3]^2)*10  # possibly anomalous electron (due to [DiffEqGPU]), intercept and cancel.
            warn_num += 1
            if warn_num < max_warn_num
                @warn "[Ensemble Simulation] Found electron (#$i in the batch) with anomalous momentum $([px,py,pz])."
            elseif warn_num == max_warn_num
                @warn "[Ensemble Simulation] Found electron (#$i in the batch) with anomalous momentum $([px,py,pz]). Similar warnings would be suppressed."
            end
            continue
        end
        phase = (traj_phase_method == :CTMC) ? (0.) : (sol.u[i][end][7])
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
            class_rates_ion[threadid] += init[8,i]
            p_inf = sqrt(2E_inf)
            a_vec = p_vec × L_vec - nucl_charge * r_vec ./ norm(r_vec)
            p_inf_vec = (p_inf/(1+p_inf^2*L2)) .* (p_inf .* (L_vec×a_vec) - a_vec)
            pxIdx = round(Int, (p_inf_vec[1]+final_p_max[1])/(final_p_max[1]/final_p_num[1]*2))
            pyIdx = round(Int, (p_inf_vec[2]+final_p_max[2])/(final_p_max[2]/final_p_num[2]*2))
            pzIdx = round(Int, (p_inf_vec[3]+final_p_max[3])/(final_p_max[3]/final_p_num[3]*2))
            if checkbounds(Bool, ion_prob_collect, pxIdx,pyIdx,pzIdx, threadid)
                if traj_phase_method == :CTMC
                    ion_prob_collect[pxIdx,pyIdx,pzIdx, threadid] += init[8,i] # ionRate
                else
                    ion_prob_collect[pxIdx,pyIdx,pzIdx, threadid] += sqrt(init[8,i])*exp(1im*init[9,i]) # sqrt(ionRate)*phaseFactor
                end
            else
                class_rates_ion_uncollected[threadid] += init[8,i]
            end
        else            # finally become rydberg.
            class_rates_ryd[threadid] += init[8,i]
            if final_ryd_collect
                n = round(Int, nucl_charge / sqrt(-2E_inf))
                l = round(Int, (sqrt(1.0+4L2)-1.0)/2)
                m = round(Int, L_vec[3])
                nIdx = n
                lIdx = l+1
                mIdx = m+final_ryd_n_max
                if checkbounds(Bool, ryd_prob_collect, nIdx,lIdx,mIdx, threadid)
                    if traj_phase_method == :CTMC
                        ryd_prob_collect[nIdx,lIdx,mIdx,threadid] += init[8,i]
                    else
                        ryd_prob_collect[nIdx,lIdx,mIdx,threadid] += sqrt(init[8,i])*exp(1im*init[9,i])
                    end
                else
                    class_rates_ryd_uncollected[threadid] += init[8,i]
                end
            else
                class_rates_ryd_uncollected[threadid] += init[8,i]
            end
        end
    end
    sum!(ion_prob_sum_temp, ion_prob_collect)
    ion_prob_final .+= ion_prob_sum_temp
    if final_ryd_collect
        sum!(ryd_prob_sum_temp, ryd_prob_collect)
        ryd_prob_final .+= ryd_prob_sum_temp
    end
    classical_rates[:ion]                += sum(class_rates_ion)
    classical_rates[:ion_uncollected]    += sum(class_rates_ion_uncollected)
    classical_rates[:ryd]                += sum(class_rates_ryd)
    classical_rates[:ryd_uncollected]    += sum(class_rates_ryd_uncollected)
end

function default_filename()
    Y,M,D = yearmonthday(now())
    h,m,s = hour(now()), minute(now()), second(now())
    return "SCSFI-$(string(Y,pad=4))$(string(M,pad=2))$(string(D,pad=2))-$(string(h,pad=2))$(string(m,pad=2))$(string(s,pad=2)).h5"
end

end
