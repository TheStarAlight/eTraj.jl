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

```@example manual_lasers
t = -600:1:600  # hide
l = Cos4Laser(peak_int=1e14, wave_len=800.0, cyc_num=10, ellip=0)   # hide
plot(t, UnitEnvelope(l).(t), color=:gray, fill=0, α=0.2)    # hide
plot!(t, -1 .* UnitEnvelope(l).(t), color=:gray, fill=0, α=0.2)     # hide
plot!(t, LaserAx(l).(t) ./ LaserA0(l), color=RGB(0/255,154/255,250/255))    # hide
plot!(xlabel=raw"$t$/a.u.", ylabel=raw"$f_\mathrm{env}$", xlim=(-600,600), ylim=(-1.05,1.05), legend=:none, size=(500,400), framestyle=:box, fmt=:svg)  # hide
```


## Cos²-envelope Laser

```@example manual_lasers
t = -600:1:600  # hide
l = Cos2Laser(peak_int=1e14, wave_len=800.0, cyc_num=10, ellip=0)   # hide
plot(t, UnitEnvelope(l).(t), color=:gray, fill=0, α=0.2)    # hide
plot!(t, -1 .* UnitEnvelope(l).(t), color=:gray, fill=0, α=0.2)     # hide
plot!(t, LaserAx(l).(t) ./ LaserA0(l), color=RGB(0/255,154/255,250/255))    # hide
plot!(xlabel=raw"$t$/a.u.", ylabel=raw"$f_\mathrm{env}$", xlim=(-600,600), ylim=(-1.05,1.05), legend=:none, size=(500,400), framestyle=:box, fmt=:svg)  # hide
```


## Gaussian-envelope Laser

```@example manual_lasers
t = -1000:1:1000  # hide
l = GaussianLaser(peak_int=1e14, wave_len=800.0, FWHM_duration=1103.1998233636713, ellip=0)   # hide
plot(t, UnitEnvelope(l).(t), color=:gray, fill=0, α=0.2)    # hide
plot!(t, -1 .* UnitEnvelope(l).(t), color=:gray, fill=0, α=0.2)     # hide
plot!(t, LaserAx(l).(t) ./ LaserA0(l), color=RGB(0/255,154/255,250/255))    # hide
plot!(xlabel=raw"$t$/a.u.", ylabel=raw"$f_\mathrm{env}$", xlim=(-1000,1000), ylim=(-1.05,1.05), legend=:none, size=(500,400), framestyle=:box, fmt=:svg)  # hide
```


## Trapezoidal-envelope Laser

```@example manual_lasers
t = -600:1:600  # hide
l = TrapezoidalLaser(peak_int=1e14, wave_len=800.0, cyc_num_turn_on=3, cyc_num_turn_off=3, cyc_num_const=4, ellip=0, t_shift=-551.6)   # hide
plot(t, UnitEnvelope(l).(t), color=:gray, fill=0, α=0.2)    # hide
plot!(t, -1 .* UnitEnvelope(l).(t), color=:gray, fill=0, α=0.2)     # hide
plot!(t, LaserAx(l).(t) ./ LaserA0(l), color=RGB(0/255,154/255,250/255))    # hide
plot!(xlabel=raw"$t$/a.u.", ylabel=raw"$f_\mathrm{env}$", xlim=(-600,600), ylim=(-1.05,1.05), legend=:none, size=(500,400), framestyle=:box, fmt=:svg)  # hide
```

