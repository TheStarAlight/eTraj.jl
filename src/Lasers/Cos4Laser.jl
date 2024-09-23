
"""
    struct Cos4Laser <: MonochromaticLaser <: Laser

Represents a monochromatic elliptically polarized laser field with Cos4-shape envelope propagating in z direction.
"""
struct Cos4Laser <: MonochromaticLaser
    peak_int;
    wave_len;
    cyc_num;
    ellip;
    azi;
    cep;
    t_shift;
end

"""
    Cos4Laser(peak_int, wave_len|ang_freq, cyc_num|duration, ellip [,azi=0.0] [,cep=0.0] [,t_shift=0.0]) <: MonochromaticLaser

Initializes a new monochromatic elliptically polarized laser field with Cos4-shape envelope.

## Parameters
- `peak_int`    : Peak intensity of the laser (numerically in **W/cm²** or a `Unitful.Quantity`).
- `wave_len`    : Wavelength of the laser (numerically in **nm** or a `Unitful.Quantity`).
- `ang_freq`    : Angular frequency of the laser (numerically in **a.u.** or a `Unitful.Quantity` of single-photon energy).
- `cyc_num`     : Number of cycles of the laser.
- `duration`    : Duration of the laser (numerically in **a.u.** or a `Unitful.Quantity`).
- `ellip`       : Ellipticity of the laser [-1≤ε≤1, 0 indicates linear polarization and ±1 indicates circular polarization].
- `azi`         : Azimuth angle of the laser's polarization's principle axis relative to x axis (numerically in radian or a `Unitful.Quantity`) *(optional, default 0)*.
- `cep`         : Carrier-Envelope-Phase of the laser (numerically in radian or a `Unitful.Quantity`) *(optional, default 0)*.
- `t_shift`     : Time shift of the laser (numerically in **a.u.** or a `Unitful.Quantity`) relative to the peak *(optional, default 0)*.

## Examples
```jldoctest
julia> l = Cos4Laser(peak_int=4e14, wave_len=800.0, cyc_num=2.0, ellip=1.0)
[MonochromaticLaser] Envelope cos⁴, peak intensity 4.0e+14 W/cm², wavelen=800 nm, 2 cycle(s), ε=1 [circularly polarized]

julia> using eTraj.Units

julia> l = Cos4Laser(peak_int=0.4PW/cm^2, ang_freq=1.5498eV, duration=5.34fs, ellip=0.0)
[MonochromaticLaser] Envelope cos⁴, peak intensity 4.0e+14 W/cm², wavelen=800.00 nm, 2.00 cycle(s), ε=0 [linearly polarized]
```
"""
function Cos4Laser(;peak_int,
                    wave_len=0, ang_freq=0,   # must specify either wave_len or ang_freq.
                    cyc_num=0,  duration=0,   # must specify either cyc_num or duration.
                    ellip, azi=0., cep=0., t_shift=0.)
    # make conversions
    (peak_int isa Quantity) && (peak_int = uconvert(W/cm^2, peak_int).val)
    (wave_len isa Quantity) && (wave_len = uconvert(nm, wave_len).val)
    (ang_freq isa Quantity) && (ang_freq = (uconvert(eV, ang_freq) |> auconvert).val)
    (duration isa Quantity) && (duration = (uconvert(fs, duration) |> auconvert).val)
    (t_shift  isa Quantity) && (t_shift  = (uconvert(fs, t_shift)  |> auconvert).val)
    (azi isa Quantity) && (azi=uconvert(u"rad",azi).val)
    (cep isa Quantity) && (cep=uconvert(u"rad",cep).val)
    # ================
    @assert wave_len>0 || ang_freq>0    "[Cos4Laser] Must specify either `wave_len` or `ang_freq`."
    @assert cyc_num>0 || duration>0     "[Cos4Laser] Must specify either `cyc_num` or `duration`."
    if wave_len>0 && ang_freq>0
        @warn "[Cos4Laser] Both `wave_len` & `ang_freq` are specified, will use `wave_len`."
    end
    if cyc_num>0 && duration>0
        @warn "[Cos4Laser] Both `cyc_num` & `duration` are specified, will use `cyc_num`."
    end
    if wave_len==0
        wave_len = 45.563352525 / ang_freq
    else
        ang_freq = 45.563352525 / wave_len
    end
    if cyc_num==0
        cyc_num = duration / (2π/ang_freq)
    end
    Cos4Laser(peak_int, wave_len, cyc_num, ellip, azi, cep, t_shift)
end

PeakInt(l::Cos4Laser) = l.peak_int
WaveLen(l::Cos4Laser) = l.wave_len
CycNum(l::Cos4Laser) = l.cyc_num
Ellipticity(l::Cos4Laser) = l.ellip
Azimuth(l::Cos4Laser) = l.azi
AngFreq(l::Cos4Laser) = 45.563352525 / l.wave_len
Period(l::Cos4Laser) = 2π / AngFreq(l)
CEP(l::Cos4Laser) = l.cep
TimeShift(l::Cos4Laser) = l.t_shift
LaserF0(l::Cos4Laser) = sqrt(l.peak_int/(1.0+l.ellip^2)/3.50944521e16)
LaserA0(l::Cos4Laser) = LaserF0(l) / AngFreq(l)
KeldyshParameter(l::Cos4Laser, Ip) = AngFreq(l) * sqrt(2Ip) / LaserF0(l)

function UnitEnvelope(l::Cos4Laser)
    local ω = AngFreq(l); local N = CycNum(l); local Δt = l.t_shift; local k = 200/Period(l)
    function (t)
        t -= Δt
        cos(ω*t/(2N))^4 * (0.5+0.5*tanh(k*(real(t)+N*π/ω))) * (0.5+0.5*tanh(-k*(real(t)-N*π/ω)))
    end
end

function LaserAx(l::Cos4Laser)
    local A0 = LaserA0(l); local ω = AngFreq(l); local N = CycNum(l); local φ = l.cep; local Δt = l.t_shift; local ε = l.ellip; local ϕ = l.azi; local k = 200/Period(l)
    return if ϕ==0
        function(t)
            t -= Δt
            A0 * cos(ω*t/(2N))^4 * (0.5+0.5*tanh(k*(real(t)+N*π/ω))) * (0.5+0.5*tanh(-k*(real(t)-N*π/ω))) * cos(ω*t+φ)
        end
    else
        function(t)
            t -= Δt
            A0 * cos(ω*t/(2N))^4 * (0.5+0.5*tanh(k*(real(t)+N*π/ω))) * (0.5+0.5*tanh(-k*(real(t)-N*π/ω))) * (cos(ω*t+φ)*cos(ϕ)+sin(ω*t+φ)*ε*sin(ϕ))
        end
    end
end
function LaserAy(l::Cos4Laser)
    local A0 = LaserA0(l); local ω = AngFreq(l); local N = CycNum(l); local φ = l.cep; local Δt = l.t_shift; local ε = l.ellip; local ϕ = l.azi; local k = 200/Period(l)
    return if ϕ==0
        function(t)
            t -= Δt
            A0 * cos(ω*t/(2N))^4 * (0.5+0.5*tanh(k*(real(t)+N*π/ω))) * (0.5+0.5*tanh(-k*(real(t)-N*π/ω))) * sin(ω*t+φ) * ε
        end
    else
        function(t)
            t -= Δt
            A0 * cos(ω*t/(2N))^4 * (0.5+0.5*tanh(k*(real(t)+N*π/ω))) * (0.5+0.5*tanh(-k*(real(t)-N*π/ω))) * (cos(ω*t+φ)*-sin(ϕ)+sin(ω*t+φ)*ε*cos(ϕ))
        end
    end
end
function LaserFx(l::Cos4Laser)
    local F0 = LaserF0(l); local ω = AngFreq(l); local N = CycNum(l); local φ = l.cep; local Δt = l.t_shift; local ε = l.ellip; local ϕ = l.azi; local k = 200/Period(l)
    return if ϕ==0
        function(t)
            t -= Δt
            F0 * cos(ω*t/(2N))^3 * (0.5+0.5*tanh(k*(real(t)+N*π/ω))) * (0.5+0.5*tanh(-k*(real(t)-N*π/ω))) * ( cos(ω*t/(2N))*sin(ω*t+φ) + 2/N*sin(ω*t/(2N))*cos(ω*t+φ))
        end
    else
        function(t)
            t -= Δt
            F0 * cos(ω*t/(2N))^3 * (0.5+0.5*tanh(k*(real(t)+N*π/ω))) * (0.5+0.5*tanh(-k*(real(t)-N*π/ω))) * ( (cos(ω*t/(2N))*sin(ω*t+φ) + 2/N*sin(ω*t/(2N))*cos(ω*t+φ))*cos(ϕ) + (-cos(ω*t/(2N))*cos(ω*t+φ) + 2/N*sin(ω*t/(2N))*sin(ω*t+φ))*ε*sin(ϕ) )
        end
    end
end
function LaserFy(l::Cos4Laser)
    local F0 = LaserF0(l); local ω = AngFreq(l); local N = CycNum(l); local φ = l.cep; local Δt = l.t_shift; local ε = l.ellip; local ϕ = l.azi; local k = 200/Period(l)
    return if ϕ==0
        function(t)
            t -= Δt
            F0 * cos(ω*t/(2N))^3 * (0.5+0.5*tanh(k*(real(t)+N*π/ω))) * (0.5+0.5*tanh(-k*(real(t)-N*π/ω))) * (-cos(ω*t/(2N))*cos(ω*t+φ) + 2/N*sin(ω*t/(2N))*sin(ω*t+φ)) * ε
        end
    else
        function(t)
            t -= Δt
            F0 * cos(ω*t/(2N))^3 * (0.5+0.5*tanh(k*(real(t)+N*π/ω))) * (0.5+0.5*tanh(-k*(real(t)-N*π/ω))) * ( (cos(ω*t/(2N))*sin(ω*t+φ) + 2/N*sin(ω*t/(2N))*cos(ω*t+φ))*-sin(ϕ) + (-cos(ω*t/(2N))*cos(ω*t+φ) + 2/N*sin(ω*t/(2N))*sin(ω*t+φ))*ε*cos(ϕ) )
        end
    end
end

function Base.show(io::IO, l::Cos4Laser)
    print(io, "[MonochromaticLaser] Envelope cos⁴, ")
    @printf(io, "peak intensity %.1e W/cm², ", l.peak_int)
    if isinteger(l.wave_len)
        @printf(io, "wavelen=%i nm, ", l.wave_len)
    else
        @printf(io, "wavelen=%.2f nm, ", l.wave_len)
    end
    if isinteger(l.cyc_num)
        @printf(io, "%i cycle(s), ", l.cyc_num)
    else
        @printf(io, "%.2f cycle(s), ", l.cyc_num)
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
    if l.t_shift != 0
        if isinteger(l.t_shift)
            @printf(io, ", peaks @ t=%i a.u.", l.t_shift)
        else
            @printf(io, ", peaks @ t=%.2f a.u.", l.t_shift)
        end
    end
    if l.cep != 0
        @printf(io, ", CEP=%.2f π", l.cep/π)
    end
    if l.azi != 0
        @printf(io, ", prin_ax_azimuth=%.2f°", l.azi/π*180)
    end
end

function Serialize(l::Cos4Laser)
    dict = OrderedDict{Symbol,Any}()
    type        = typeof(l)
    peak_int    = l.peak_int
    wave_len    = l.wave_len
    cyc_num     = l.cyc_num
    ellip       = l.ellip
    azi         = l.azi
    cep         = l.cep
    t_shift     = l.t_shift
    @pack! dict = (type, peak_int, wave_len, cyc_num, ellip, azi, cep, t_shift)
    return dict
end
