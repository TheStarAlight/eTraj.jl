# Lasers

*This section provides information of available lasers in the library.*

In this section we list available lasers implemented in the [`Lasers`](@ref) module of the library.

```@docs
Lasers
```

```@contents
Pages = ["manual2_lasers.md"]
Depth = 3
```


```@setup manual_lasers
using SemiclassicalSFI
using SemiclassicalSFI.Lasers
using Plots
pyplot()
```

## Cos⁴-envelope Laser

A [`Lasers.Cos4Laser`](@ref) has a cos⁴-shaped-envelope:
```math
f_{\mathrm{env}}(t) =
    \begin{cases}
    \cos^4{\left[ \omega (t-t_0)/2N \right]},    & -NT/2 \leq t-t_0 \leq NT/2, \\
    0,                                          & \mathrm{otherwise,} \\
    \end{cases}
```
where ``\omega`` is the angular frequency of the laser field, ``N`` is the cycle number, ``T`` the period, and ``t_0`` denotes the peak time (corresponding to the `time_shift` in the constructor method).

For the detailed usage, cf. the documentation of [`Lasers.Cos4Laser`](@ref).

```@docs
Lasers.Cos4Laser
```

The following shows an example of [`Lasers.Cos4Laser`](@ref) and its envelope shape.

```@repl manual_lasers
Cos4Laser(peak_int=1e14, wave_len=800.0, cyc_num=10, ellip=0)
```

```@example manual_lasers
t = -600:1:600  # hide
l = Cos4Laser(peak_int=1e14, wave_len=800.0, cyc_num=10, ellip=0)   # hide
plot(t, UnitEnvelope(l).(t), color=:gray, fill=0, α=0.2)    # hide
plot!(t, -1 .* UnitEnvelope(l).(t), color=:gray, fill=0, α=0.2)     # hide
plot!(t, LaserAx(l).(t) ./ LaserA0(l), color=RGB(0/255,154/255,250/255))    # hide
plot!(xlabel=raw"$t$/a.u.", ylabel=raw"$f_\mathrm{env}$", xlim=(-600,600), ylim=(-1.05,1.05), legend=:none, size=(500,400), framestyle=:box, fmt=:svg)  # hide
```


## Cos²-envelope Laser

A [`Lasers.Cos2Laser`](@ref) has a cos²-shaped-envelope, similiar with the `Cos4Laser`:
```math
f_{\mathrm{env}}(t) =
    \begin{cases}
    \cos^2{\left[ \omega (t-t_0)/2N \right]},   & -NT/2 \leq t-t_0 \leq NT/2, \\
    0,                                          & \mathrm{otherwise.} \\
    \end{cases}
```

For the detailed usage, cf. the documentation of [`Lasers.Cos2Laser`](@ref).

```@docs
Lasers.Cos2Laser
```

The following shows an example of [`Lasers.Cos2Laser`](@ref) and its envelope shape.

```@repl manual_lasers
Cos2Laser(peak_int=1e14, wave_len=800.0, cyc_num=10, ellip=0)
```

```@example manual_lasers
t = -600:1:600  # hide
l = Cos2Laser(peak_int=1e14, wave_len=800.0, cyc_num=10, ellip=0)   # hide
plot(t, UnitEnvelope(l).(t), color=:gray, fill=0, α=0.2)    # hide
plot!(t, -1 .* UnitEnvelope(l).(t), color=:gray, fill=0, α=0.2)     # hide
plot!(t, LaserAx(l).(t) ./ LaserA0(l), color=RGB(0/255,154/255,250/255))    # hide
plot!(xlabel=raw"$t$/a.u.", ylabel=raw"$f_\mathrm{env}$", xlim=(-600,600), ylim=(-1.05,1.05), legend=:none, size=(500,400), framestyle=:box, fmt=:svg)  # hide
```


## Gaussian-envelope Laser

A [`Lasers.GaussianLaser`](@ref) has a Gaussian-shaped-envelope:
```math
f_{\mathrm{env}}(t) = \exp{\left[ -(t-t_0)^2/\sigma^2 \right]} = \exp{\left[ -8\ln{2}(t-t_0)^2/\tau_{\mathrm{FWHM}}^2 \right]},
```
where ``\sigma`` is the temporal width of the laser (relating to the `spread_duration` in the constructor method) and ``\tau_{\mathrm{FWHM}}=2\sqrt{2\ln{2}}\sigma`` denotes the laser's temporal FWHM (Full-Width at Half Maximum).

For the detailed usage, cf. the documentation of [`Lasers.GaussianLaser`](@ref).

```@docs
Lasers.GaussianLaser
```

An example of [`Lasers.GaussianLaser`](@ref) and its envelope shape are shown as follows.

```@repl manual_lasers
GaussianLaser(peak_int=1e14, wave_len=800.0, FWHM_duration=1103.2, ellip=0)
```

```@example manual_lasers
t = -1000:1:1000  # hide
l = GaussianLaser(peak_int=1e14, wave_len=800.0, FWHM_duration=1103.2, ellip=0)   # hide
plot(t, UnitEnvelope(l).(t), color=:gray, fill=0, α=0.2)    # hide
plot!(t, -1 .* UnitEnvelope(l).(t), color=:gray, fill=0, α=0.2)     # hide
plot!(t, LaserAx(l).(t) ./ LaserA0(l), color=RGB(0/255,154/255,250/255))    # hide
plot!(xlabel=raw"$t$/a.u.", ylabel=raw"$f_\mathrm{env}$", xlim=(-1000,1000), ylim=(-1.05,1.05), legend=:none, size=(500,400), framestyle=:box, fmt=:svg)  # hide
```


## Trapezoidal-envelope Laser

A [`Lasers.TrapezoidalLaser`](@ref) has a trapezoidal-shaped-envelope:
```math
f_{\mathrm{env}}(t) =
    \begin{cases}
        (t-t_0) / N_{\mathrm{on}}T,
                & 0 < t-t_0 \leq N_{\mathrm{on}}T, \\
        1,
                & N_{\mathrm{on}}T < t-t_0 \leq (N_{\mathrm{on}}+N_{\mathrm{const}})T, \\
        1-[t-t_0-(N_{\mathrm{on}}+N_{\mathrm{const}})T] / N_{\mathrm{off}}T,
                & (N_{\mathrm{on}}+N_{\mathrm{const}})T < t-t_0 \leq (N_{\mathrm{on}}+N_{\mathrm{const}}+N_{\mathrm{off}})T, \\
        0,
                & \mathrm{otherwise.}
    \end{cases}
```
where ``N_{\mathrm{on}}, N_{\mathrm{const}}, N_{\mathrm{off}}`` are cycle numbers during the turn-on, constant, and turn-off stages. *Note that for the `TrapezoidalLaser`, the ``t_0`` denotes the time of rise instead of time of peak, in contrast to the previous lasers.*

For the detailed usage, cf. the documentation of [`Lasers.TrapezoidalLaser`](@ref).

```@docs
Lasers.TrapezoidalLaser
```

In the following is an example of [`Lasers.TrapezoidalLaser`](@ref) and its envelope shape.

```@repl manual_lasers
l = TrapezoidalLaser(
        peak_int=1e14, wave_len=800.0,
        cyc_num_turn_on=3, cyc_num_turn_off=3, cyc_num_const=4,
        ellip=0, t_shift=-551.6)
```

```@example manual_lasers
t = -600:1:600  # hide
l = TrapezoidalLaser(peak_int=1e14, wave_len=800.0, cyc_num_turn_on=3, cyc_num_turn_off=3, cyc_num_const=4, ellip=0, t_shift=-551.6)   # hide
plot(t, UnitEnvelope(l).(t), color=:gray, fill=0, α=0.2)    # hide
plot!(t, -1 .* UnitEnvelope(l).(t), color=:gray, fill=0, α=0.2)     # hide
plot!(t, LaserAx(l).(t) ./ LaserA0(l), color=RGB(0/255,154/255,250/255))    # hide
plot!(xlabel=raw"$t$/a.u.", ylabel=raw"$f_\mathrm{env}$", xlim=(-600,600), ylim=(-1.05,1.05), legend=:none, size=(500,400), framestyle=:box, fmt=:svg)  # hide
```

