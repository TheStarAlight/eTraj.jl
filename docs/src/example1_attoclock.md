# Example: Attoclock & Initial Condition Methods

This example is adapted from [[*JPB* **54**, 144001 (2021)](https://doi.org/10.1088/1361-6455/ac0d3e)].

In the present example of attoclock experiment, an ultra-short pulse of circular polarization is employed to explore the ultrafast attosecond dynamics such as the tunneling time delay.
In this example we would perform simulation of an attoclock experiment using different initial condition methods, namely ADK, SFA-SPANE and SFA-SPA, which reveals the influence of non-adiabatic effects to the attoclock signal.
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

The momentum spectra are shown in the figure below.
Due to the exponential dependence of ionization rate on the field strength, the PMD exhibits a crescent shape near the peak of the negative vector potential ``-\AA(t)``.
In an adiabatic tunneling scenario, which corresponds to the ADK initial condition, the trace of ``-\AA(t)`` is expected to match the median of the crescent shape.
While for non-adiabatic tunneling, the distribution of the initial transverse momentum ``\kkt`` at the tunnel exit centers at a nonzero value, which results in expansion of the crescent shape and the enhancement of overall ionization probability, as can be seen in the figure.
Moreover, the PMD obtained with the SFA-SPANE and SFA-SPA initial conditions shows similar shape and total ionization probability, which demonstrate the SFA-SPANE's advantage of preserving the non-adiabatic effects with much less computational cost than the SFA-SPA.

![fig:example_2cycs_CP](assets/figure_2cycs_CP.png)
