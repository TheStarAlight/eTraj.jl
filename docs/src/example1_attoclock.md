# Example: Attoclock & Initial Condition Methods

!!! note "Guide on Generating Plots"
    We attached scripts for generating plots for each example in the `examples/` directory.
    This script requires additional installation of the `Plots.jl` and `PyPlot.jl` packages.
    To install the dependencies, first, follow the instructions [here](@ref config_py) to configure the `PyCall.jl` package and install the `matplotlib` python package.
    Then, install `Plots.jl` and `PyPlot.jl` packages by running the following commands:
    ```julia
    using Pkg
    Pkg.add("Plots")
    Pkg.add("PyPlot")
    ```
    After installing the dependencies, you can run the plot scripts either in a Julia REPL or directly from the command line, see the ["Running Scripts" section](@ref running-scripts).

This example is adapted from [[*JPB* **54**, 144001 (2021)](https://doi.org/10.1088/1361-6455/ac0d3e)].

The attoclock experiment employs an ultra-short circularly polarized pulse to investigate ultrafast attosecond dynamics, particularly tunneling time delay phenomena.
Here, we simulate an attoclock experiment using three initial condition methods—ADK, SFA-SPANE, and SFA-SPA—to examine how non-adiabatic effects influence the attoclock signal.
We choose CTMC as the phase method because the quantum interference effect is not significant in this example.
In our theoretical framework, the simulation schemes are named after "ADK-CTMC", "SFA-SPANE-CTMC" and "SFA-SPA-CTMC", respectively.

```julia
# examples/test_2cycs_CP.jl

using eTraj
using eTraj.Targets, eTraj.Lasers, eTraj.Units

l = Cos4Laser(peak_int=0.4PW/cm^2, wave_len=800.0nm, cyc_num=2, ellip=1.0)
t = get_atom("H")

for init_cond in [:ADK, :SPANE, :SPA]
    perform_traj_simulation(
        init_cond_method    = init_cond,
        traj_phase_method   = :CTMC,
        laser               = l,
        target              = t,
        dimension           = 2,            # 2D simulation, x-y plane only
        sample_t_intv       = (-100,100),   # equivalent to `(-2.42fs, 2.42fs)`
        sample_t_num        = 20000,        # will sample 20000 equidistant time points between -100 and 100 a.u.
        traj_t_final        = 120,          # the traj end at 120 a.u., equivalent to `2.90fs`
        final_p_max         = (2.5,2.5),    # the momentum spec collection grid's border (-2.5 to +2.5 a.u.)
        final_p_num         = (500,500),    # the momentum spec collection grid's size (500x500)
        ss_kd_max           = 2.0,
        ss_kd_num           = 10000,        # will sample 10000 equidistant k⟂ points between -2 to +2 a.u.
        output_path         = "$(init_cond)-CTMC_4e14_800nm_cos4_2cyc_CP.jld2"
    )
end
```

The PMD results are presented in the figure below.
Owing to the exponential dependence of the ionization rate on the field strength, the PMD exhibits a crescent-shaped structure near the peak of the negative vector potential ``-\AA(t)``.
In the adiabatic tunneling regime, corresponding to the ADK initial condition, the trajectory of ``-\AA(t)`` is expected to align with the median of the crescent-shaped structure.
In contrast, for non-adiabatic tunneling, the distribution of the initial transverse momentum ``\kkt`` at the tunnel exit is centered at a nonzero value. This leads to an expansion of the crescent-shaped structure and an enhancement of the overall ionization probability, as evident in the figure.
Furthermore, the PMDs obtained using the SFA-SPANE and SFA-SPA initial conditions exhibit similar shapes and total ionization probabilities. This demonstrates the advantage of the SFA-SPANE approach in preserving non-adiabatic effects while significantly reducing computational costs compared to the SFA-SPA method.

![fig:example_2cycs_CP](assets/figure_2cycs_CP.png)
