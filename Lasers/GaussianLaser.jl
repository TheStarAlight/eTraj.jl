
"Represents a monochromatic elliptically polarized laser field with Gaussian-shape envelope propagating in z direction."
struct GaussianLaser <: MonochromaticLaser
    "Peak intensity of the laser field (in W/cm^2)."
    peak_int;
    "Wavelength of the laser field (in NANOMETER)."
    wave_len;
    "Temporal width (converting to cycle numbers) of the laser field, namely σ."
    spread_cyc_num;
    "Ellipticity of the laser field."
    ellip;
    "Azimuth angle of the laser's polarization's principle axis relative to x axis (in radians)."
    azi;
    "Carrier-Envelope-Phase (CEP) of the laser field."
    cep;
    "Time shift of the laser relative to the peak (in a.u.)."
    t_shift;
    """
    Constructs a new monochromatic elliptically polarized laser field with Gaussian-shape envelope.
    # Parameters
    - `peakInt`     : Peak intensity of the laser field (in W/cm²).
    - `WaveLen`     : Wavelength of the laser field (in nm).
    - `spreadCycNum`: Temporal width (converting to cycle numbers) of the laser field, namely σ.
    - `ellip`       : Ellipticity of the laser field [-1≤e≤1, 0 indicates linear polarization and ±1 indicates circular polarization].
    - `azi`         : Azimuth angle of the laser's polarization's principle axis relative to x axis (in radians) (optional, default 0).
    - `cep`         : Carrier-Envelope-Phase of the laser field (optional, default 0).
    - `t_shift`     : Time shift of the laser (in a.u.) relative to the peak (optional, default 0).
    """
    function GaussianLaser(peakInt, waveLen, spreadCycNum, ellip, azi=0., cep=0., t_shift=0.)
        @assert peakInt>0       "[GaussianLaser] Peak intensity must be positive."
        @assert waveLen>0       "[GaussianLaser] Wavelength must be positive."
        @assert spreadCycNum>0  "[GaussianLaser] Cycle number must be positive."
        @assert -1≤ellip≤1      "[GaussianLaser] Ellipticity must be in [-1,1]."
        new(peakInt,waveLen,spreadCycNum,ellip,azi,cep,t_shift)
    end
    """
    Constructs a new monochromatic elliptically polarized laser field with Gaussian-shape envelope.
    # Parameters
    - `peakInt`         : Peak intensity of the laser field (in W/cm²).
    - `WaveLen`         : Wave length of the laser field (in nm). Must specify either `waveLen` or `angFreq`.
    - `angFreq`         : Angular frequency of the laser field (in a.u.). Must specify either `waveLen` or `angFreq`.
    - `spreadCycNum`    : Temporal width (converting to cycle numbers) of the laser field, namely σ. Must specify one in `spreadCycNum`, `spreadDuration` and `FWHM_duration`.
    - `spreadDuration`  : Temporal width of the laser field (in a.u.). Must specify one in `spreadCycNum`, `spreadDuration` and `FWHM_duration`.
    - `FWHM_duration`   : Temporal FWHM(Full Width at Half Maxima) of the laser field (in a.u.). Must specify one in `spreadCycNum`, `spreadDuration` and `FWHM_duration`.
    - `ellip`           : Ellipticity of the laser field [-1≤e≤1, 0 indicates linear polarization and ±1 indicates circular polarization].
    - `azi`             : Azimuth angle of the laser's polarization's principle axis relative to x axis (in radians) (optional, default 0).
    - `cep`             : Carrier-Envelope-Phase of the laser field (optional, default 0).
    - `t_shift`         : Time shift of the laser (in a.u.) relative to the peak (optional, default 0).
    """
    function GaussianLaser(;peakInt,
                            waveLen=-1, angFreq=-1,     # must specify either waveLen or angFreq.
                            spreadCycNum=-1, spreadDuration=-1, FWHM_duration=-1,   # must specify one in spreadCycNum, spreadDuration and FWHM_duration.
                            ellip, azi=0., cep=0., t_shift=0.)
        @assert waveLen>0 || angFreq>0                                  "[GaussianLaser] Must specify either waveLen or angFreq."
        @assert spreadCycNum>0 || spreadDuration>0 || FWHM_duration>0   "[GaussianLaser] Must specify one in spreadCycNum, spreadDuration and FWHM_duration."
        if waveLen>0 && angFreq>0
            @warn "[GaussianLaser] Both waveLen & angFreq are specified, will use waveLen."
        end
        if 1*(spreadCycNum>0) + 1*(spreadDuration>0) + 1*(FWHM_duration>0) > 1
            @warn "[GaussianLaser] More than one in spreadCycNum & spreadDuration & FWHM_duration are specified, will use the first positive value."
        end
        if waveLen==-1
            waveLen = 45.563352525 / angFreq
        end
        if spreadCycNum==-1
            if spreadDuration==-1
                spreadCycNum = FWHM_duration / (2*sqrt(2*log(2))) / (2π/angFreq)
            else
                spreadCycNum = spreadDuration / (2π/angFreq)
            end
        end
        GaussianLaser(peakInt,waveLen,spreadCycNum,ellip,azi,cep,t_shift)
    end
end
"Gets the peak intensity of the laser field (in W/cm²)."
PeakInt(l::GaussianLaser) = l.peak_int
"Gets the wave length of the laser field (in nm)."
WaveLen(l::GaussianLaser) = l.wave_len
"Gets the temporal width (converting to cycle numbers) of the laser field, namely σ."
SpreadCycNum(l::GaussianLaser) = l.spread_cyc_num
"Gets the temporal width (in a.u.) of the laser field, namely σ."
SpreadDuration(l::GaussianLaser) = l.spread_cyc_num * Period(l)
"Gets the temporal FWHM(Full Width at Half Maxima) of the laser field (in a.u.)."
FWHM_Duration(l::GaussianLaser) = l.spread_cyc_num * Period(l) * (2*sqrt(2*log(2)))
"Gets the ellipticity of the laser field."
Ellipticity(l::GaussianLaser) = l.ellip
"Gets the azimuth angle of the laser's polarization's principle axis relative to x axis (in radians)."
Azimuth(l::GaussianLaser) = l.azi
"Gets the angular frequency (ω) of the laser field (in a.u.)."
AngFreq(l::GaussianLaser) = 45.563352525 / l.wave_len
"Gets the period of the laser field (in a.u.)."
Period(l::GaussianLaser) = 2π / AngFreq(l)
"Gets the time shift relative to the peak (in a.u.)."
TimeShift(l::GaussianLaser) = l.t_shift
"Gets the peak electric field intensity of the laser field (in a.u.)."
LaserF0(l::GaussianLaser) = sqrt(l.peak_int/(1.0+l.ellip^2)/3.50944521e16)
"Gets the peak vector potential intensity of the laser field (in a.u.)."
LaserA0(l::GaussianLaser) = LaserF0(l) / AngFreq(l)

"Gets the time-dependent x component of the vector potential under dipole approximation."
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
"Gets the time-dependent y component of the vector potential under dipole approximation."
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
"Gets the time-dependent x component of the electric field strength under dipole approximation."
function LaserFx(l::GaussianLaser)
    local A0 = LaserA0(l); local ω = AngFreq(l); local σ = SpreadDuration(l); local φ = l.cep; local Δt = l.t_shift; local ε = l.ellip; local ϕ = l.azi;
    return if ϕ==0
        function(t)
            t -= Δt
            A0 * exp(-t^2/2/σ^2) * (-t/σ^2*cos(ω*t+φ) + ω*sin(ω*t+φ))
        end
    else
        function(t)
            t -= Δt
            A0 * exp(-t^2/2/σ^2) * (-(t/σ^2*cos(ω*t+φ) - ω*sin(ω*t+φ))*cos(ϕ) + (t/σ^2*sin(ω*t+φ)-ω*cos(ω*t+φ))*ε*sin(ϕ))
        end
    end
end
"Gets the time-dependent y component of the electric field strength under dipole approximation."
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
            A0 * exp(-t^2/2/σ^2) * (-(t/σ^2*cos(ω*t+φ) - ω*sin(ω*t+φ))*-sin(ϕ) + (t/σ^2*sin(ω*t+φ)-ω*cos(ω*t+φ))*ε*cos(ϕ))
        end
    end
end

"Prints the information about the laser."
Base.show(io::IO, l::GaussianLaser) = print(io,"[MonochromaticLaser] Envelope Gaussian, Wavelength=$(l.waveLen) nm, Temporal width $(l.spreadCycNum) cycle(s) [FWHM $(FWHM_Duration(l)*24.19e-3) fs], e=$(l.ellip)"
                                                * (l.ellip==0 ? " [Linearly polarized]" : "") * (abs(l.ellip)==1 ? " [Circularly polarized]" : "")
                                                * ", PrincipleAxisAzimuth=$(l.azi/π*180)°" * (l.t_shift==0 ? "" : ", Peaks at t₀=$(l.t_shift) a.u.") * (l.cep==0 ? "" : ", CEP=$(l.cep)"))