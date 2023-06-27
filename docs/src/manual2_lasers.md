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

```@meta
CurrentModule = SemiclassicalSFI.Lasers
```

```@setup manual_lasers
using SemiclassicalSFI
using SemiclassicalSFI.Lasers
using Plots
pyplot()
```


## Basic Properties of Monochromatic Lasers

A monochromatic laser is composed of the carrier wave ``\cos{(\omega t+\phi)}`` and the envelope ``f_{\mathrm{env}}(t)``.
Given the amplitude of the vector potential ``A_0``, the time-dependent vector potential of the laser,
which we assume to propagate in ``z`` direction and have ``x`` axis as the principle axis of polarization, reads
```math
\bm{A}(t) =  A_0 f_{\mathrm{env}}(t) \cos{(\omega t+\phi)} \left( \bm{e}_x + \varepsilon \bm{e}_y \right),
```
where ``\varepsilon`` is the ellipticity.

### List of Available Properties

Currently the monochromatic lasers implemented in the library include
[`Cos4Laser`](@ref), [`Cos2Laser`](@ref), [`GaussianLaser`](@ref) and [`TrapezoidalLaser`](@ref).

The available properties of the laser fields are listed below.
To obtain a property of the laser field, invoke the property as a method and pass the laser object as an argument. The following shows an example.

```@repl manual_lasers
l = Cos4Laser(peak_int=1e14, wave_len=800.0, cyc_num=10, ellip=0)
LaserA0(l)
Ax = LaserAx(l)
Ax(0.0)
```

|               |[`Cos4Laser`](@ref) | [`Cos2Laser`](@ref) | [`GaussianLaser`](@ref) | [`TrapezoidalLaser`](@ref) |
|:--------------|:-:|:-:|:-:|:-:|
|`PeakInt`      | ✔ | ✔ | ✔ | ✔ |
|`WaveLen`      | ✔ | ✔ | ✔ | ✔ |
|`CycNum`       | ✔ | ✔ |   |   |
|`SpreadCycNum` |   |   | ✔ |   |
|`SpreadDuration`|  |   | ✔ |   |
|`FWHM_Duration`|   |   | ✔ |   |
|`CycNumTotal`  |   |   |   | ✔ |
|`CycNumTurnOn` |   |   |   | ✔ |
|`CycNumTurnOff`|   |   |   | ✔ |
|`CycNumConst`  |   |   |   | ✔ |
|`Ellipticity`  | ✔ | ✔ | ✔ | ✔ |
|`Azimuth`      | ✔ | ✔ | ✔ | ✔ |
|`AngFreq`      | ✔ | ✔ | ✔ | ✔ |
|`Period`       | ✔ | ✔ | ✔ | ✔ |
|`CEP`          | ✔ | ✔ | ✔ | ✔ |
|`TimeShift`    | ✔ | ✔ | ✔ | ✔ |
|`LaserF0`      | ✔ | ✔ | ✔ | ✔ |
|`LaserA0`      | ✔ | ✔ | ✔ | ✔ |
|`LaserFx`      | ✔ | ✔ | ✔ | ✔ |
|`LaserFy`      | ✔ | ✔ | ✔ | ✔ |
|`LaserAx`      | ✔ | ✔ | ✔ | ✔ |
|`LaserAy`      | ✔ | ✔ | ✔ | ✔ |

### Ellipticity

The ellipticity ``\varepsilon`` defines the polarization type of the laser field.
For special cases, `0` indicates linear polarization and `±1` indicates circular polarization. The electric field rotates clockwise for positive ellipticities and counter-clockwise for negative ones.

```@setup manual_lasers
l = Cos4Laser(peak_int=1e14, wave_len=800.0, cyc_num=6, ellip=0)
F0 = LaserF0(l)
Fx = LaserFx(l)
Fy = LaserFy(l)
t = -400:1:400
Fxt = Fx.(t)/F0
Fyt = Fy.(t)/F0
color = RGB(0/255,154/255,250/255)
color_trans = RGBA(0/255,154/255,250/255,0.3)
f1=plot(t,Fxt,Fyt, label="ε=0.0", color=color)
plot!(repeat([-400],length(t)),Fxt,Fyt, label=:none, color=color_trans)
plot!(t,Fxt,repeat([-1.0],length(t)), label=:none, color=color_trans)
plot!(t,repeat([1.0],length(t)),Fyt, label=:none, color=color_trans)
plot!(xlim=(-400,400),ylim=(-1,1),zlim=(-1,1))
plot!(xlabel=raw"$t$/a.u.", ylabel=raw"$F_x/F_0$", zlabel=raw"$F_y/F_0$")

l = Cos4Laser(peak_int=1e14, wave_len=800.0, cyc_num=6, ellip=0.5)
F0 = LaserF0(l)
Fx = LaserFx(l)
Fy = LaserFy(l)
t = -400:1:400
Fxt = Fx.(t)/F0
Fyt = Fy.(t)/F0
color = RGB(255/255,127/255,24/255)
color_trans = RGBA(255/255,127/255,24/255,0.3)
f2=plot(t,Fxt,Fyt, label="ε=0.5", color=color)
plot!(repeat([-400],length(t)),Fxt,Fyt, label=:none, color=color_trans)
plot!(t,Fxt,repeat([-1.0],length(t)), label=:none, color=color_trans)
plot!(t,repeat([1.0],length(t)),Fyt, label=:none, color=color_trans)
plot!(xlim=(-400,400),ylim=(-1,1),zlim=(-1,1))
plot!(xlabel=raw"$t$/a.u.", ylabel=raw"$F_x/F_0$", zlabel=raw"$F_y/F_0$")

l = Cos4Laser(peak_int=1e14, wave_len=800.0, cyc_num=6, ellip=1.0)
F0 = LaserF0(l)
Fx = LaserFx(l)
Fy = LaserFy(l)
t = -400:1:400
Fxt = Fx.(t)/F0
Fyt = Fy.(t)/F0
color = RGB(44/255,160/255,44/255)
color_trans = RGBA(44/255,160/255,44/255,0.3)
f3=plot(t,Fxt,Fyt, label="ε=1.0", color=color)
plot!(repeat([-400],length(t)),Fxt,Fyt, label=:none, color=color_trans)
plot!(t,Fxt,repeat([-1.0],length(t)), label=:none, color=color_trans)
plot!(t,repeat([1.0],length(t)),Fyt, label=:none, color=color_trans)
plot!(xlim=(-400,400),ylim=(-1,1),zlim=(-1,1))
plot!(xlabel=raw"$t$/a.u.", ylabel=raw"$F_x/F_0$", zlabel=raw"$F_y/F_0$")

l = Cos4Laser(peak_int=1e14, wave_len=800.0, cyc_num=6, ellip=-1.0)
F0 = LaserF0(l)
Fx = LaserFx(l)
Fy = LaserFy(l)
t = -400:1:400
Fxt = Fx.(t)/F0
Fyt = Fy.(t)/F0
color = RGB(214/255,39/255,40/255)
color_trans = RGBA(214/255,39/255,40/255,0.3)
f4=plot(t,Fxt,Fyt, label="ε=-1.0", color=color)
plot!(repeat([-400],length(t)),Fxt,Fyt, label=:none, color=color_trans)
plot!(t,Fxt,repeat([-1.0],length(t)), label=:none, color=color_trans)
plot!(t,repeat([1.0],length(t)),Fyt, label=:none, color=color_trans)
plot!(xlim=(-400,400),ylim=(-1,1),zlim=(-1,1))
plot!(xlabel=raw"$t$/a.u.", ylabel=raw"$F_x/F_0$", zlabel=raw"$F_y/F_0$")

f=plot(f1,f2,f3,f4, layout=(2,2), size=(1000,800), fmt=:svg)
savefig("./manual2_lasers_ellip.svg")
```
![manual2_lasers_ellip.svg](manual2_lasers_ellip.svg)


### Azimuth of Principle Axis

The azimuth angle ``\varphi`` of the principle axis defines a clockwise rotation of the laser field in the polarization plane.

```@setup manual_lasers
l = Cos4Laser(peak_int=1e14, wave_len=800.0, cyc_num=6, ellip=0, azi=0π)
F0 = LaserF0(l)
Fx = LaserFx(l)
Fy = LaserFy(l)
t = -400:1:400
Fxt = Fx.(t)/F0
Fyt = Fy.(t)/F0
color = RGB(0/255,154/255,250/255)
color_trans = RGBA(0/255,154/255,250/255,0.3)
f1=plot(t,Fxt,Fyt, label="azi=0", color=color)
plot!(repeat([-400],length(t)),Fxt,Fyt, label=:none, color=color_trans)
plot!(t,Fxt,repeat([-1.0],length(t)), label=:none, color=color_trans)
plot!(t,repeat([1.0],length(t)),Fyt, label=:none, color=color_trans)
plot!(xlim=(-400,400),ylim=(-1,1),zlim=(-1,1))
plot!(xlabel=raw"$t$/a.u.", ylabel=raw"$F_x/F_0$", zlabel=raw"$F_y/F_0$")

l = Cos4Laser(peak_int=1e14, wave_len=800.0, cyc_num=6, ellip=0, azi=π/6)
F0 = LaserF0(l)
Fx = LaserFx(l)
Fy = LaserFy(l)
t = -400:1:400
Fxt = Fx.(t)/F0
Fyt = Fy.(t)/F0
color = RGB(255/255,127/255,24/255)
color_trans = RGBA(255/255,127/255,24/255,0.3)
f2=plot(t,Fxt,Fyt, label="azi=π/6 (30°)", color=color)
plot!(repeat([-400],length(t)),Fxt,Fyt, label=:none, color=color_trans)
plot!(t,Fxt,repeat([-1.0],length(t)), label=:none, color=color_trans)
plot!(t,repeat([1.0],length(t)),Fyt, label=:none, color=color_trans)
plot!(xlim=(-400,400),ylim=(-1,1),zlim=(-1,1))
plot!(xlabel=raw"$t$/a.u.", ylabel=raw"$F_x/F_0$", zlabel=raw"$F_y/F_0$")

l = Cos4Laser(peak_int=1e14, wave_len=800.0, cyc_num=6, ellip=0, azi=π/3)
F0 = LaserF0(l)
Fx = LaserFx(l)
Fy = LaserFy(l)
t = -400:1:400
Fxt = Fx.(t)/F0
Fyt = Fy.(t)/F0
color = RGB(44/255,160/255,44/255)
color_trans = RGBA(44/255,160/255,44/255,0.3)
f3=plot(t,Fxt,Fyt, label="azi=π/3 (60°)", color=color)
plot!(repeat([-400],length(t)),Fxt,Fyt, label=:none, color=color_trans)
plot!(t,Fxt,repeat([-1.0],length(t)), label=:none, color=color_trans)
plot!(t,repeat([1.0],length(t)),Fyt, label=:none, color=color_trans)
plot!(xlim=(-400,400),ylim=(-1,1),zlim=(-1,1))
plot!(xlabel=raw"$t$/a.u.", ylabel=raw"$F_x/F_0$", zlabel=raw"$F_y/F_0$")

l = Cos4Laser(peak_int=1e14, wave_len=800.0, cyc_num=6, ellip=0, azi=π/2)
F0 = LaserF0(l)
Fx = LaserFx(l)
Fy = LaserFy(l)
t = -400:1:400
Fxt = Fx.(t)/F0
Fyt = Fy.(t)/F0
color = RGB(214/255,39/255,40/255)
color_trans = RGBA(214/255,39/255,40/255,0.3)
f4=plot(t,Fxt,Fyt, label="azi=π/2 (90°)", color=color)
plot!(repeat([-400],length(t)),Fxt,Fyt, label=:none, color=color_trans)
plot!(t,Fxt,repeat([-1.0],length(t)), label=:none, color=color_trans)
plot!(t,repeat([1.0],length(t)),Fyt, label=:none, color=color_trans)
plot!(xlim=(-400,400),ylim=(-1,1),zlim=(-1,1))
plot!(xlabel=raw"$t$/a.u.", ylabel=raw"$F_x/F_0$", zlabel=raw"$F_y/F_0$")

f=plot(f1,f2,f3,f4, layout=(2,2), size=(1000,800), fmt=:svg)
savefig("./manual2_lasers_azimuth.svg")
```
![manual2_lasers_azimuth.svg](manual2_lasers_azimuth.svg)


### Carrier-Envelope-Phase (CEP)

The carrier-envelope-phase (CEP) ``\phi`` is the difference between the optical phase of the carrier wave and the envelope position. For few-cycle laser pulses, the influence of the CEP to the laser-matter interaction becomes significant.

```@setup manual_lasers
t = -400:1:400
l = Cos4Laser(peak_int=1e14, wave_len=800.0, cyc_num=6, ellip=0, cep=0)
plot(t, UnitEnvelope(l).(t), color=:gray, fill=0, α=0.1, label=:none)
plot!(t, -1 .* UnitEnvelope(l).(t), color=:gray, fill=0, α=0.1, label=:none)

l = Cos4Laser(peak_int=1e14, wave_len=800.0, cyc_num=6, ellip=0, cep=0)
plot!(t, LaserAx(l).(t) ./ LaserA0(l), label="CEP=0", color=RGB(64/255,183/255,173/255))
l = Cos4Laser(peak_int=1e14, wave_len=800.0, cyc_num=6, ellip=0, cep=π/3)
plot!(t, LaserAx(l).(t) ./ LaserA0(l), label="CEP=π/3", color=RGB(52/255,143/255,167/255))
l = Cos4Laser(peak_int=1e14, wave_len=800.0, cyc_num=6, ellip=0, cep=2π/3)
plot!(t, LaserAx(l).(t) ./ LaserA0(l), label="CEP=2π/3", color=RGB(55/255,101/255,158/255))
l = Cos4Laser(peak_int=1e14, wave_len=800.0, cyc_num=6, ellip=0, cep=π)
plot!(t, LaserAx(l).(t) ./ LaserA0(l), label="CEP=π", color=RGB(65/255,61/255,163/255))

f=plot!(xlabel=raw"$t$/a.u.", ylabel=raw"$f_\mathrm{env}$", xlim=(-300,300), ylim=(-1.05,1.05), size=(500,400), framestyle=:box, fmt=:svg)
savefig("./manual2_lasers_cep.svg")
```
![manual2_lasers_cep.svg](manual2_lasers_cep.svg)


## Cos⁴-envelope Laser

A [`Cos4Laser`](@ref) has a cos⁴-shaped-envelope:
```math
f_{\mathrm{env}}(t) =
    \begin{cases}
    \cos^4{\left[ \omega (t-t_0)/2N \right]},   & -NT/2 \leq t-t_0 \leq NT/2, \\
    0,                                          & \mathrm{otherwise,} \\
    \end{cases}
```
where ``\omega`` is the angular frequency of the laser field, ``N`` is the cycle number, ``T`` the period, and ``t_0`` denotes the peak time (corresponding to the `time_shift` in the constructor method).

For the detailed usage, cf. the documentation of [`Cos4Laser`](@ref).

```@docs
Lasers.Cos4Laser
```

The following shows an example of [`Cos4Laser`](@ref) and its envelope shape.

```@repl manual_lasers
Cos4Laser(peak_int=1e14, wave_len=800.0, cyc_num=10, ellip=0)
```

```@setup manual_lasers
t = -600:1:600
l = Cos4Laser(peak_int=1e14, wave_len=800.0, cyc_num=10, ellip=0)
plot(t, UnitEnvelope(l).(t), color=:gray, fill=0, α=0.2)
plot!(t, -1 .* UnitEnvelope(l).(t), color=:gray, fill=0, α=0.2)
plot!(t, LaserAx(l).(t) ./ LaserA0(l), color=RGB(0/255,154/255,250/255))
plot!(xlabel=raw"$t$/a.u.", ylabel=raw"$f_\mathrm{env}$", xlim=(-600,600), ylim=(-1.05,1.05), legend=:none, size=(500,400), framestyle=:box, fmt=:svg)
savefig("./manual2_lasers_env_cos4.svg")
```
![manual2_lasers_env_cos4.svg](manual2_lasers_env_cos4.svg)

## Cos²-envelope Laser

A [`Cos2Laser`](@ref) has a cos²-shaped-envelope, similiar with the `Cos4Laser`:
```math
f_{\mathrm{env}}(t) =
    \begin{cases}
    \cos^2{\left[ \omega (t-t_0)/2N \right]},   & -NT/2 \leq t-t_0 \leq NT/2, \\
    0,                                          & \mathrm{otherwise.} \\
    \end{cases}
```

For the detailed usage, cf. the documentation of [`Cos2Laser`](@ref).

```@docs
Lasers.Cos2Laser
```

The following shows an example of [`Cos2Laser`](@ref) and its envelope shape.

```@repl manual_lasers
Cos2Laser(peak_int=1e14, wave_len=800.0, cyc_num=10, ellip=0)
```

```@setup manual_lasers
t = -600:1:600
l = Cos2Laser(peak_int=1e14, wave_len=800.0, cyc_num=10, ellip=0)
plot(t, UnitEnvelope(l).(t), color=:gray, fill=0, α=0.2)
plot!(t, -1 .* UnitEnvelope(l).(t), color=:gray, fill=0, α=0.2)
plot!(t, LaserAx(l).(t) ./ LaserA0(l), color=RGB(0/255,154/255,250/255))
plot!(xlabel=raw"$t$/a.u.", ylabel=raw"$f_\mathrm{env}$", xlim=(-600,600), ylim=(-1.05,1.05), legend=:none, size=(500,400), framestyle=:box, fmt=:svg)
savefig("./manual2_lasers_env_cos2.svg")
```
![manual2_lasers_env_cos2.svg](manual2_lasers_env_cos2.svg)


## Gaussian-envelope Laser

A [`GaussianLaser`](@ref) has a Gaussian-shaped-envelope:
```math
f_{\mathrm{env}}(t) = \exp{\left[ -(t-t_0)^2/\sigma^2 \right]} = \exp{\left[ -8\ln{2}(t-t_0)^2/\tau_{\mathrm{FWHM}}^2 \right]},
```
where ``\sigma`` is the temporal width of the laser (relating to the `spread_duration` in the constructor method) and ``\tau_{\mathrm{FWHM}}=2\sqrt{2\ln{2}}\sigma`` denotes the laser's temporal FWHM (Full-Width at Half Maximum).

For the detailed usage, cf. the documentation of [`GaussianLaser`](@ref).

```@docs
Lasers.GaussianLaser
```

An example of [`GaussianLaser`](@ref) and its envelope shape are shown as follows.

```@repl manual_lasers
GaussianLaser(peak_int=1e14, wave_len=800.0, FWHM_duration=1103.2, ellip=0)
```

```@setup manual_lasers
t = -1000:1:1000
l = GaussianLaser(peak_int=1e14, wave_len=800.0, FWHM_duration=1103.2, ellip=0)
plot(t, UnitEnvelope(l).(t), color=:gray, fill=0, α=0.2)
plot!(t, -1 .* UnitEnvelope(l).(t), color=:gray, fill=0, α=0.2)
plot!(t, LaserAx(l).(t) ./ LaserA0(l), color=RGB(0/255,154/255,250/255))
plot!(xlabel=raw"$t$/a.u.", ylabel=raw"$f_\mathrm{env}$", xlim=(-1000,1000), ylim=(-1.05,1.05), legend=:none, size=(500,400), framestyle=:box, fmt=:svg)
savefig("./manual2_lasers_env_gaussian.svg")
```
![manual2_lasers_env_gaussian.svg](manual2_lasers_env_gaussian.svg)


## Trapezoidal-envelope Laser

A [`TrapezoidalLaser`](@ref) has a trapezoidal-shaped-envelope:
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
where ``N_{\mathrm{on}}, N_{\mathrm{const}}, N_{\mathrm{off}}`` are cycle numbers during the turn-on, constant, and turn-off stages. Note that for the `TrapezoidalLaser`, the ``t_0`` denotes the *time of rise* instead of time of peak, in contrast to the previous lasers.

For the detailed usage, cf. the documentation of [`TrapezoidalLaser`](@ref).

```@docs
Lasers.TrapezoidalLaser
```

In the following is an example of [`TrapezoidalLaser`](@ref) and its envelope shape.

```@repl manual_lasers
l = TrapezoidalLaser(
        peak_int=1e14, wave_len=800.0,
        cyc_num_turn_on=3, cyc_num_turn_off=3, cyc_num_const=4,
        ellip=0, t_shift=-551.6)
```

```@setup manual_lasers
t = -600:1:600
l = TrapezoidalLaser(peak_int=1e14, wave_len=800.0, cyc_num_turn_on=3, cyc_num_turn_off=3, cyc_num_const=4, ellip=0, t_shift=-551.6)
plot(t, UnitEnvelope(l).(t), color=:gray, fill=0, α=0.2)
plot!(t, -1 .* UnitEnvelope(l).(t), color=:gray, fill=0, α=0.2)
plot!(t, LaserAx(l).(t) ./ LaserA0(l), color=RGB(0/255,154/255,250/255))
plot!(xlabel=raw"$t$/a.u.", ylabel=raw"$f_\mathrm{env}$", xlim=(-600,600), ylim=(-1.05,1.05), legend=:none, size=(500,400), framestyle=:box, fmt=:svg)
savefig("./manual2_lasers_env_trapezoidal.svg")
```
![manual2_lasers_env_trapezoidal.svg](manual2_lasers_env_trapezoidal.svg)

