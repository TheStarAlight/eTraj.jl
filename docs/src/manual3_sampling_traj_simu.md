# [Manual: Electron Sampling and Trajectory Simulation](@id traj_simu)

The `ElectronSamplers` module provides means of generating initial electron samples using different initial condition methods.
The `ElectronSampler` is an abstract supertype, with `ADKSampler`, `SPANESampler`, `SPASampler` and `WFATSampler` being its subtypes.
When the user starts a trajectory simulation job by invoking [`eTraj.perform_traj_simulation`](@ref), the method would further call the `init_sampler` method, which would assign the corresponding type of `ElectronSampler` in the background.
These invocations execute in the background thus the `ElectronSamplers` module is kept internal (i.e., not exported for public invocation).

The [`eTraj.perform_traj_simulation`](@ref) method serves as a public entrance to performing a trajectory simulation.
The method would automatically detect number of available threads (specified by passing command-line arguments `-t <thread_num>` when starting julia) and run the trajectory simulation in parallel.

Here we brief on the working procedure of the method [`eTraj.perform_traj_simulation`](@ref):

- First, the `eTraj.perform_traj_simulation` method initializes an `eTraj.TrajectorySimulationJob`, which stores the essential parameters. The electron sampler is assigned according to `init_cond_method`.
- Then, it repeatedly invokes the `eTraj.launch_and_collect!` method, where in each invocation a batch of electrons which tunnel at the same time ``\tr`` but have different transversal momenta ``\kkt`` is launched and collected using the corresponding simulation scheme.
  - The sampling behavior is controlled by the `sample_monte_carlo` parameter, if `sample_monte_carlo` is `true`, the ``\tr`` of the electron batches and the ``\kkt`` in each batch would be randomly sampled inside the given intervals specified by `sample_t_intv`, `mc_kd_max` and `mc_kz_max`,
  otherwise, the initial conditions of the electrons would be sampled in an equidistant manner, with the intervals controlled by `sample_t_intv`, `ss_kd_max` and `ss_kz_max`.
- After generating a batch of initial electrons, the `eTraj.launch_and_collect!` method simulates the electrons' classical trajectories (together with the phase), and then collects the electrons' final momenta on the grid. Trajectory simulations are carried out using the [`OrdinaryDiffEq.jl`](https://github.com/SciML/OrdinaryDiffEq.jl) package in the [`DifferentialEquations.jl`](https://github.com/SciML/DifferentialEquations.jl) ecosystem.
  - The grid's size and spacing depends on the `final_p_max` and `final_p_num` parameters.
- Finally, the `eTraj.perform_traj_simulation` method generates the output file which contains the photoelectron momentum distribution (PMD) and other necessary information in a Julia Data Format (JLD2) file or an HDF5 file.

The output file is a Julia JLD2 file, which is compatible with the HDF5 data format, and can be opened by the [`JLD2.jl`](https://github.com/JuliaIO/JLD2.jl) package:
```julia-repl
julia> using JLD2

julia> file = jldopen("ADK-CTMC_4e14_800nm_cos4_2cyc_CP.jld2")
JLDFile .../ADK-CTMC_4e14_800nm_cos4_2cyc_CP.jld2 (read-only)
    ├─ info
    ├─ params_text              # parameters stored in YAML format
    ├─ params                   # parameters stored in a Julia `Dict`
    ├─ px                       # coordinates of the `momentum_spec` on x-axis
    ├─ py                       # coordinates of the `momentum_spec` on y-axis
    ├─ momentum_spec            # PMD data stored in a Julia `Array`
    ├─ ion_prob                 # total ionization probability
    ├─ ion_prob_uncollected     # ionization probability of discarded electrons
    └─ num_effective_traj       # total number of effective trajectories

julia> file["params_text"] |> print  # equivalent to print(file["params_text"])
init_cond_method: ADK
laser:
    type: Cos4Laser
    peak_int: 4.0e14
    wave_len: 800.0
    ⋮
target:
    type: HydrogenLikeAtom
    Ip: 0.5
    nucl_charge: 1
    ⋮
sample_t_intv: (-100, 100)
sample_t_num: 20000
sample_monte_carlo: false
⋮
```

The detailed documentation of the [`eTraj.perform_traj_simulation`](@ref) is shown below:

```@docs
eTraj.perform_traj_simulation
```

