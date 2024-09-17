
"""
    struct TrapezoidalLaser <: MonochromaticLaser <: Laser

Represents a monochromatic elliptically polarized laser field with Trapezoidal-shape envelope propagating in z direction.
"""
struct TrapezoidalLaser <: MonochromaticLaser
    peak_int;
    wave_len;
    cyc_num_turn_on;
    cyc_num_turn_off;
    cyc_num_const;
    ellip;
    azi;
    cep;
    t_turn_on;
end

"""
    TrapezoidalLaser(peak_int, wave_len|ang_freq, cyc_num_turn_on, cyc_num_turn_off, cyc_num_const, ellip [,azi=0.0] [,cep=0.0] [,t_turn_on=0.0]) <: MonochromaticLaser

Initializes a new monochromatic elliptically polarized laser field with trapezoidal-shape envelope.

## Parameters
- `peak_int`        : Peak intensity of the laser field (numerically in **W/cm²** or a `Unitful.Quantity`).
- `wave_len`        : Wavelength of the laser field (numerically in **nm** or a `Unitful.Quantity`).
- `ang_freq`        : Angular frequency of the laser field (numerically in **a.u.** or a `Unitful.Quantity` of single-photon energy).
- `cyc_num_turn_on` : Number of cycles of the laser field in the turn-on.
- `cyc_num_turn_off`: Number of cycles of the laser field in the turn-off.
- `cyc_num_const`   : Number of cycles of the laser field in the constant-intensity.
- `ellip`           : Ellipticity of the laser field [-1≤ε≤1, 0 indicates linear polarization and ±1 indicates circular polarization].
- `azi`             : Azimuth angle of the laser's polarization's principle axis relative to x axis (numerically in radian or a `Unitful.Quantity`) *(optional, default 0)*.
- `cep`             : Carrier-Envelope-Phase of the laser field (numerically in radian or a `Unitful.Quantity`) *(optional, default 0)*.
- `t_turn_on`       : Time shift of the laser (numerically in **a.u.** or a `Unitful.Quantity`) relative to the beginning of TURN-ON *(optional, default 0)*.

## Examples
```jldoctest
julia> l = TrapezoidalLaser(peak_int=4e14, wave_len=800.0, cyc_num_turn_on=2, cyc_num_turn_off=2, cyc_num_const=6, ellip=1.0)
[MonochromaticLaser] Envelope Trapezoidal, peak intensity 4.0e+14 W/cm², wavelen=800 nm, turn_on/const/turn_off 2/6/2 cycle(s), ε=1 [circularly polarized]

julia> using SemiclassicalSFI.Units

julia> TrapezoidalLaser(peak_int=4e14W/cm^2, wave_len=800.0nm, cyc_num_turn_on=2, cyc_num_turn_off=2, cyc_num_const=6, ellip=0.0, t_turn_on=-10.0fs)
[MonochromaticLaser] Envelope Trapezoidal, peak intensity 4.0e+14 W/cm², wavelen=800 nm, turn_on/const/turn_off 2/6/2 cycle(s), ε=0 [linearly polarized], rises @ t=-413.41 a.u.
```
"""
function TrapezoidalLaser(; peak_int,
                            wave_len=0, ang_freq=0,     # must specify either wave_len or ang_freq.
                            cyc_num_turn_on, cyc_num_turn_off, cyc_num_const,
                            ellip, azi=0., cep=0., t_turn_on=0.)
    # make conversions
    (peak_int isa Quantity) && (peak_int = uconvert(W/cm^2, peak_int).val)
    (wave_len isa Quantity) && (wave_len = uconvert(nm, wave_len).val)
    (ang_freq isa Quantity) && (ang_freq = (uconvert(eV, ang_freq) |> auconvert).val)
    (t_turn_on isa Quantity) && (t_turn_on = (uconvert(fs, t_turn_on) |> auconvert).val)
    (azi isa Quantity) && (azi=uconvert(u"rad",azi).val)
    (cep isa Quantity) && (cep=uconvert(u"rad",cep).val)
    # ================
    @assert wave_len>0 || ang_freq>0    "[TrapezoidalLaser] Must specify either `wave_len` or `ang_freq`."
    @assert cyc_num_turn_on>0 && cyc_num_turn_off>0 && cyc_num_const≥0  "[TrapezoidalLaser] Cycle numbers must be positive."
    if wave_len>0 && ang_freq>0
        @warn "[TrapezoidalLaser] Both `wave_len` & `ang_freq` are specified, will use `wave_len`."
    end
    if wave_len==0
        wave_len = 45.563352525 / ang_freq
    end
    TrapezoidalLaser(peak_int, wave_len, cyc_num_turn_on, cyc_num_turn_off, cyc_num_const, ellip, azi, cep, t_turn_on)
end

PeakInt(l::TrapezoidalLaser) = l.peak_int
WaveLen(l::TrapezoidalLaser) = l.wave_len
CycNumTotal(l::TrapezoidalLaser) = l.cyc_num_turn_on + l.cyc_num_turn_off + l.cyc_num_const
CycNumTurnOn(l::TrapezoidalLaser) = l.cyc_num_turn_on
CycNumTurnOff(l::TrapezoidalLaser) = l.cyc_num_turn_off
CycNumConst(l::TrapezoidalLaser) = l.cyc_num_const
Ellipticity(l::TrapezoidalLaser) = l.ellip
Azimuth(l::TrapezoidalLaser) = l.azi
AngFreq(l::TrapezoidalLaser) = 45.563352525 / l.wave_len
Period(l::TrapezoidalLaser) = 2π / AngFreq(l)
CEP(l::TrapezoidalLaser) = l.cep
TimeTurnOn(l::TrapezoidalLaser) = l.t_turn_on
LaserF0(l::TrapezoidalLaser) = sqrt(l.peak_int/(1.0+l.ellip^2)/3.50944521e16)
LaserA0(l::TrapezoidalLaser) = LaserF0(l) / AngFreq(l)
KeldyshParameter(l::TrapezoidalLaser, Ip) = AngFreq(l) * sqrt(2Ip) / LaserF0(l)

function UnitEnvelope(l::TrapezoidalLaser)
    local T = Period(l); local Δt = TimeShift(l);
    local N_on = CycNumTurnOn(l); local N_const = CycNumConst(l); local N_off = CycNumTurnOff(l);
    local t_on = N_on*T; local t_const = N_const*T; local t_off = N_off*T;
    function (t)
        t -= Δt
        ( (0<t<t_on)*(t/t_on) + (t_on≤t≤(t_on+t_const))*1.0 + ((t_on+t_const)<t<(t_on+t_const+t_off))*(1-(t-t_on-t_const)/t_off) )
    end
end

function LaserAx(l::TrapezoidalLaser)
    local A0 = LaserA0(l); local ω = AngFreq(l); local T = Period(l); local Δt = TimeShift(l);
    local N_on = CycNumTurnOn(l); local N_const = CycNumConst(l); local N_off = CycNumTurnOff(l);
    local t_on = N_on*T; local t_const = N_const*T; local t_off = N_off*T;
    local φ = l.cep; local ε = l.ellip; local ϕ = l.azi;
    return if ϕ==0
        function(t)  # this function might be accessed by GPU, thus the conditional expressions should be avoided.
            t -= Δt
            A0 * ( (0<t<t_on)*(t/t_on) + (t_on≤t≤(t_on+t_const))*1.0 + ((t_on+t_const)<t<(t_on+t_const+t_off))*(1-(t-t_on-t_const)/t_off) ) * cos(ω*t+φ)
        end
    else
        function(t)
            t -= Δt
            A0 * ( (0<t<t_on)*(t/t_on) + (t_on≤t≤(t_on+t_const))*1.0 + ((t_on+t_const)<t<(t_on+t_const+t_off))*(1-(t-t_on-t_const)/t_off) )*(cos(ω*t+φ)*cos(ϕ) + sin(ω*t+φ)*ε*sin(ϕ))
        end
    end
end
function LaserAy(l::TrapezoidalLaser)
    local A0 = LaserA0(l); local ω = AngFreq(l); local T = Period(l); local Δt = TimeShift(l);
    local N_on = CycNumTurnOn(l); local N_const = CycNumConst(l); local N_off = CycNumTurnOff(l);
    local t_on = N_on*T; local t_const = N_const*T; local t_off = N_off*T;
    local φ = l.cep; local ε = l.ellip; local ϕ = l.azi;
    return if ϕ==0
        function(t)
            t -= Δt
            A0 * ( (0<t<t_on)*(t/t_on) + (t_on≤t≤(t_on+t_const))*1.0 + ((t_on+t_const)<t<(t_on+t_const+t_off))*(1-(t-t_on-t_const)/t_off) ) * sin(ω*t+φ) * ε
        end
    else
        function(t)
            t -= Δt
            A0 * ( (0<t<t_on)*(t/t_on) + (t_on≤t≤(t_on+t_const))*1.0 + ((t_on+t_const)<t<(t_on+t_const+t_off))*(1-(t-t_on-t_const)/t_off) )*(cos(ω*t+φ)*-sin(ϕ) + sin(ω*t+φ)*ε*cos(ϕ))
        end
    end
end
function LaserFx(l::TrapezoidalLaser)
    local A0 = LaserA0(l); local ω = AngFreq(l); local T = Period(l); local Δt = TimeShift(l);
    local N_on = CycNumTurnOn(l); local N_const = CycNumConst(l); local N_off = CycNumTurnOff(l);
    local t_on = N_on*T; local t_const = N_const*T; local t_off = N_off*T;
    local φ = l.cep; local ε = l.ellip; local ϕ = l.azi;
    return if ϕ==0
        function(t)
            t -= Δt
            A0 * ( ((0<t<t_on)*(t/t_on)+(t_on≤t≤(t_on+t_const))*1.0+((t_on+t_const)<t<(t_on+t_const+t_off))*(1-(t-t_on-t_const)/t_off))*sin(ω*t+φ)*ω - ((0<t<t_on)*(1/t_on)+((t_on+t_const)<t<(t_on+t_const+t_off))*(-1/t_off))*cos(ω*t+φ))
        end
    else
        function(t)
            t -= Δt
            A0 * ( (((0<t<t_on)*(t/t_on)+(t_on≤t≤(t_on+t_const))*1.0+((t_on+t_const)<t<(t_on+t_const+t_off))*(1-(t-t_on-t_const)/t_off))*sin(ω*t+φ)*ω + ((0<t<t_on)*(1/t_on)+((t_on+t_const)<t<(t_on+t_const+t_off))*(-1/t_off))*cos(ω*t+φ))*cos(ϕ) + (((0<t<t_on)*(t/t_on)+(t_on≤t≤(t_on+t_const))*1.0+((t_on+t_const)<t<(t_on+t_const+t_off))*(1-(t-t_on-t_const)/t_off))*cos(ω*t+φ)*ω + ((0<t<t_on)*(1/t_on)+((t_on+t_const)<t<(t_on+t_const+t_off))*(-1/t_off))*sin(ω*t+φ))*-ε*sin(ϕ) )
        end
    end
end
function LaserFy(l::TrapezoidalLaser)
    local A0 = LaserA0(l); local ω = AngFreq(l); local T = Period(l); local Δt = TimeShift(l);
    local N_on = CycNumTurnOn(l); local N_const = CycNumConst(l); local N_off = CycNumTurnOff(l);
    local t_on = N_on*T; local t_const = N_const*T; local t_off = N_off*T;
    local φ = l.cep; local ε = l.ellip; local ϕ = l.azi;
    return if ϕ==0
        function(t)
            t -= Δt
            -A0 * ε * ( ((0<t<t_on)*(t/t_on)+(t_on≤t≤(t_on+t_const))*1.0+((t_on+t_const)<t<(t_on+t_const+t_off))*(1-(t-t_on-t_const)/t_off))*cos(ω*t+φ)*ω + ((0<t<t_on)*(1/t_on)+((t_on+t_const)<t<(t_on+t_const+t_off))*(-1/t_off))*sin(ω*t+φ))
        end
    else
        function(t)
            t -= Δt
            A0 * ( (((0<t<t_on)*(t/t_on)+(t_on≤t≤(t_on+t_const))*1.0+((t_on+t_const)<t<(t_on+t_const+t_off))*(1-(t-t_on-t_const)/t_off))*sin(ω*t+φ)*ω + ((0<t<t_on)*(1/t_on)+((t_on+t_const)<t<(t_on+t_const+t_off))*(-1/t_off))*cos(ω*t+φ))*-sin(ϕ) + (((0<t<t_on)*(t/t_on)+(t_on≤t≤(t_on+t_const))*1.0+((t_on+t_const)<t<(t_on+t_const+t_off))*(1-(t-t_on-t_const)/t_off))*cos(ω*t+φ)*ω + ((0<t<t_on)*(1/t_on)+((t_on+t_const)<t<(t_on+t_const+t_off))*(-1/t_off))*sin(ω*t+φ))*-ε*cos(ϕ) )
        end
    end
end

function Base.show(io::IO, l::TrapezoidalLaser)
    print(io, "[MonochromaticLaser] Envelope Trapezoidal, ")
    @printf(io, "peak intensity %.1e W/cm², ", l.peak_int)
    if isinteger(l.wave_len)
        @printf(io, "wavelen=%i nm, ", l.wave_len)
    else
        @printf(io, "wavelen=%.2f nm, ", l.wave_len)
    end
    if isinteger(l.cyc_num_turn_on) && isinteger(l.cyc_num_const) && isinteger(l.cyc_num_turn_off)
        @printf(io, "turn_on/const/turn_off %i/%i/%i cycle(s), ", l.cyc_num_turn_on, l.cyc_num_const, l.cyc_num_turn_off)
    else
        @printf(io, "turn_on/const/turn_off %.2f/%.2f/%.2f cycle(s), ", l.cyc_num_turn_on, l.cyc_num_const, l.cyc_num_turn_off)
    end
    if isinteger(l.ellip)
        @printf(io, "ε=%i", l.ellip)
    else
        @printf(io, "ε=%.2f", l.ellip)
    end
    if l.ellip == 0
        print(io, " [linearly polarized]")
    elseif abs(l.ellip) == 1
        print(io, " [circularly polarized]")
    end
    if l.t_turn_on != 0
        if isinteger(l.t_turn_on)
            @printf(io, ", rises @ t=%i a.u.", l.t_turn_on)
        else
            @printf(io, ", rises @ t=%.2f a.u.", l.t_turn_on)
        end
    end
    if l.cep != 0
        @printf(io, ", CEP=%.2f π", l.cep/π)
    end
    if l.azi != 0
        @printf(io, ", prin_ax_azimuth=%.2f°", l.azi/π*180)
    end
end

function Serialize(l::TrapezoidalLaser)
    dict = OrderedDict{Symbol,Any}()
    type                = typeof(l)
    peak_int            = l.peak_int
    wave_len            = l.wave_len
    cyc_num_turn_on     = l.cyc_num_turn_on
    cyc_num_const       = l.cyc_num_const
    cyc_num_turn_off    = l.cyc_num_turn_off
    ellip               = l.ellip
    azi                 = l.azi
    cep                 = l.cep
    t_shift             = l.t_turn_on
    @pack! dict = (type, peak_int, wave_len, cyc_num_turn_on, cyc_num_const, cyc_num_turn_off, ellip, azi, cep, t_shift)
    return dict
end
