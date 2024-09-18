
"""
    struct GaussianLaser <: MonochromaticLaser <: Laser

Represents a monochromatic elliptically polarized laser field with Gaussian-shape envelope propagating in z direction.
"""
struct GaussianLaser <: MonochromaticLaser
    peak_int;
    wave_len;
    spread_cyc_num;
    ellip;
    azi;
    cep;
    t_shift;
end

"""
    GaussianLaser(peak_int, wave_len|ang_freq, spread_cyc_num|spread_duration|FWHM_duration, ellip [,azi=0.0] [,cep=0.0] [,t_shift=0.0]) <: MonochromaticLaser

Initializes a new monochromatic elliptically polarized laser field with Gaussian-shape envelope.

## Parameters
- `peak_int`        : Peak intensity of the laser field (numerically in **W/cm²** or a `Unitful.Quantity`).
- `wave_len`        : Wavelength of the laser field (numerically in **nm** or a `Unitful.Quantity`).
- `ang_freq`        : Angular frequency of the laser field (numerically in **a.u.** or a `Unitful.Quantity` of single-photon energy).
- `spread_cyc_num`  : Half-temporal width (converted to **cycle numbers**) of the laser field, namely σ.
- `spread_duration` : Half-temporal width of the laser field (numerically in **a.u.** or a `Unitful.Quantity`).
- `FWHM_duration`   : Temporal FWHM (Full Width at Half Maxima) of the laser field (numerically in **a.u.** or a `Unitful.Quantity`).
- `ellip`           : Ellipticity of the laser field [-1≤ε≤1, 0 indicates linear polarization and ±1 indicates circular polarization].
- `azi`             : Azimuth angle of the laser's polarization's principle axis relative to x axis (numerically in radian or a `Unitful.Quantity`) *(optional, default 0)*.
- `cep`             : Carrier-Envelope-Phase of the laser field (numerically in radian or a `Unitful.Quantity`) *(optional, default 0)*.
- `t_shift`         : Time shift of the laser (numerically in **a.u.** or a `Unitful.Quantity`) relative to the peak *(optional, default 0)*.

## Examples
```jldoctest
julia> l = GaussianLaser(peak_int=4e14, wave_len=800.0, spread_cyc_num=2.0, ellip=1.0)
[MonochromaticLaser] Envelope Gaussian, peak intensity 4.0e+14 W/cm², wavelen=800 nm, temporal width 4 cycle(s) [FWHM 12.57 fs], ε=1 [circularly polarized]

julia> using SemiclassicalSFI.Units

julia> l = GaussianLaser(peak_int=0.4PW/cm^2, ang_freq=1.5498eV, FWHM_duration=12.57fs, ellip=0.0)
[MonochromaticLaser] Envelope Gaussian, peak intensity 4.0e+14 W/cm², wavelen=800.00 nm, temporal width 4.00 cycle(s) [FWHM 12.57 fs], ε=0 [linearly polarized]
```
"""
function GaussianLaser(;peak_int,
                        wave_len=0, ang_freq=0,   # must specify either wave_len or ang_freq.
                        spread_cyc_num=0, spread_duration=0, FWHM_duration=0,   # must specify one in spread_cyc_num, spread_duration and FWHM_duration.
                        ellip, azi=0., cep=0., t_shift=0.)
    # make conversions
    (peak_int isa Quantity) && (peak_int = uconvert(W/cm^2, peak_int).val)
    (wave_len isa Quantity) && (wave_len = uconvert(nm, wave_len).val)
    (ang_freq isa Quantity) && (ang_freq = (uconvert(eV, ang_freq) |> auconvert).val)
    (spread_duration isa Quantity) && (spread_duration = (uconvert(fs, spread_duration) |> auconvert).val)
    (FWHM_duration isa Quantity) && (FWHM_duration = (uconvert(fs, FWHM_duration) |> auconvert).val)
    (t_shift  isa Quantity) && (t_shift  = (uconvert(fs, t_shift) |> auconvert).val)
    (azi isa Quantity) && (azi=uconvert(u"rad",azi).val)
    (cep isa Quantity) && (cep=uconvert(u"rad",cep).val)
    # ================
    @assert wave_len>0 || ang_freq>0                                    "[GaussianLaser] Must specify either `wave_len` or `ang_freq`."
    @assert spread_cyc_num>0 || spread_duration>0 || FWHM_duration>0    "[GaussianLaser] Must specify one in `spread_cyc_num`, `spread_duration` and `FWHM_duration`."
    if wave_len>0 && ang_freq>0
        @warn "[GaussianLaser] Both `wave_len` & `ang_freq` are specified, will use `wave_len`."
    end
    if 1*(spread_cyc_num>0) + 1*(spread_duration>0) + 1*(FWHM_duration>0) > 1
        @warn "[GaussianLaser] More than one in `spread_cyc_num` & `spread_duration` & `FWHM_duration` are specified, will use the first positive value."
    end
    if wave_len==0
        wave_len = 45.563352525 / ang_freq
    else
        ang_freq = 45.563352525 / wave_len
    end
    if spread_cyc_num==0
        if spread_duration==0
            spread_cyc_num = FWHM_duration / (2*sqrt(2*log(2))) / (2π/ang_freq)
        else
            spread_cyc_num = spread_duration / (2π/ang_freq)
        end
    end
    GaussianLaser(peak_int, wave_len, spread_cyc_num, ellip, azi, cep, t_shift)
end

PeakInt(l::GaussianLaser) = l.peak_int
WaveLen(l::GaussianLaser) = l.wave_len
SpreadCycNum(l::GaussianLaser) = l.spread_cyc_num
SpreadDuration(l::GaussianLaser) = l.spread_cyc_num * Period(l)
FWHM_Duration(l::GaussianLaser) = l.spread_cyc_num * Period(l) * (2*sqrt(2*log(2)))
Ellipticity(l::GaussianLaser) = l.ellip
Azimuth(l::GaussianLaser) = l.azi
AngFreq(l::GaussianLaser) = 45.563352525 / l.wave_len
Period(l::GaussianLaser) = 2π / AngFreq(l)
CEP(l::GaussianLaser) = l.cep
TimeShift(l::GaussianLaser) = l.t_shift
LaserF0(l::GaussianLaser) = sqrt(l.peak_int/(1.0+l.ellip^2)/3.50944521e16)
LaserA0(l::GaussianLaser) = LaserF0(l) / AngFreq(l)
KeldyshParameter(l::GaussianLaser, Ip) = AngFreq(l) * sqrt(2Ip) / LaserF0(l)

function UnitEnvelope(l::GaussianLaser)
    local ω = AngFreq(l); local σ = SpreadDuration(l); local Δt = l.t_shift;
    function (t)
        t -= Δt
        exp(-t^2/2/σ^2)
    end
end

function LaserAx(l::GaussianLaser)
    local A0 = LaserA0(l); local ω = AngFreq(l); local σ = SpreadDuration(l); local φ = l.cep; local Δt = l.t_shift; local ε = l.ellip; local ϕ = l.azi;
    return if ϕ==0
        function(t)
            t -= Δt
            A0 * exp(-t^2/2/σ^2) * cos(ω*t+φ)
        end
    else
        function(t)
            t -= Δt
            A0 * exp(-t^2/2/σ^2) * (cos(ω*t+φ)*cos(ϕ)+sin(ω*t+φ)*ε*sin(ϕ))
        end
    end
end
function LaserAy(l::GaussianLaser)
    local A0 = LaserA0(l); local ω = AngFreq(l); local σ = SpreadDuration(l); local φ = l.cep; local Δt = l.t_shift; local ε = l.ellip; local ϕ = l.azi;
    return if ϕ==0
        function(t)
            t -= Δt
            A0 * exp(-t^2/2/σ^2) * sin(ω*t+φ) * ε
        end
    else
        function(t)
            t -= Δt
            A0 * exp(-t^2/2/σ^2) * (cos(ω*t+φ)*-sin(ϕ)+sin(ω*t+φ)*ε*cos(ϕ))
        end
    end
end
function LaserFx(l::GaussianLaser)
    local A0 = LaserA0(l); local ω = AngFreq(l); local σ = SpreadDuration(l); local φ = l.cep; local Δt = l.t_shift; local ε = l.ellip; local ϕ = l.azi;
    return if ϕ==0
        function(t)
            t -= Δt
            A0 * exp(-t^2/2/σ^2) * (t/σ^2*cos(ω*t+φ) + ω*sin(ω*t+φ))
        end
    else
        function(t)
            t -= Δt
            A0 * exp(-t^2/2/σ^2) * ((t/σ^2*cos(ω*t+φ) + ω*sin(ω*t+φ))*cos(ϕ) + (t/σ^2*sin(ω*t+φ) - ω*cos(ω*t+φ))*ε*sin(ϕ))
        end
    end
end
function LaserFy(l::GaussianLaser)
    local A0 = LaserA0(l); local ω = AngFreq(l); local σ = SpreadDuration(l); local φ = l.cep; local Δt = l.t_shift; local ε = l.ellip; local ϕ = l.azi;
    return if ϕ==0
        function(t)
            t -= Δt
            A0 * exp(-t^2/2/σ^2) * (t/σ^2*sin(ω*t+φ) - ω*cos(ω*t+φ)) * ε
        end
    else
        function(t)
            t -= Δt
            A0 * exp(-t^2/2/σ^2) * ((t/σ^2*cos(ω*t+φ) + ω*sin(ω*t+φ))*-sin(ϕ) + (t/σ^2*sin(ω*t+φ) - ω*cos(ω*t+φ))*ε*cos(ϕ))
        end
    end
end

function Base.show(io::IO, l::GaussianLaser)
    print(io, "[MonochromaticLaser] Envelope Gaussian, ")
    @printf(io, "peak intensity %.1e W/cm², ", l.peak_int)
    if isinteger(l.wave_len)
        @printf(io, "wavelen=%i nm, ", l.wave_len)
    else
        @printf(io, "wavelen=%.2f nm, ", l.wave_len)
    end
    if isinteger(l.spread_cyc_num)
        @printf(io, "temporal width %i cycle(s) [FWHM %.2f fs], ", 2*l.spread_cyc_num, FWHM_Duration(l)*24.19e-3)
    else
        @printf(io, "temporal width %.2f cycle(s) [FWHM %.2f fs], ", 2*l.spread_cyc_num, FWHM_Duration(l)*24.19e-3)
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

function Serialize(l::GaussianLaser)
    dict = OrderedDict{Symbol,Any}()
    type            = typeof(l)
    peak_int        = l.peak_int
    wave_len        = l.wave_len
    spread_cyc_num  = l.spread_cyc_num
    ellip           = l.ellip
    azi             = l.azi
    cep             = l.cep
    t_shift         = l.t_shift
    @pack! dict = (type, peak_int, wave_len, spread_cyc_num, ellip, azi, cep, t_shift)
    return dict
end
