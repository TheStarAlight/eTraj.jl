
"Represents a monochromatic elliptically polarized laser field with Trapezoidal-shape envelope propagating in z direction."
struct TrapezoidalLaser <: MonochromaticLaser
    "Peak intensity of the laser field (in W/cm^2)."
    peak_int;
    "Wavelength of the laser field (in NANOMETER)."
    wave_len;
    "Cycle number of the laser field in the turn-on."
    cyc_num_turn_on;
    "Cycle number of the laser field in the turn-off."
    cyc_num_turn_off;
    "Cycle number of the laser field in the constant-intensity."
    cyc_num_const;
    "Ellipticity of the laser field."
    ellip;
    "Azimuth angle of the laser's polarization's principle axis relative to x axis (in radians)."
    azi;
    "Carrier-Envelope-Phase (CEP) of the laser field."
    cep;
    "Time shift of the laser relative to the beginning of TURN-ON (in a.u.)."
    t_shift;
    """
    Constructs a new monochromatic elliptically polarized laser field with Trapezoidal-shape envelope.
    # Parameters
    - `peakInt`         : Peak intensity of the laser field (in W/cm²).
    - `WaveLen`         : Wavelength of the laser field (in nm).
    - `cycNumTurnOn`    : Number of cycles of the laser field in the turn-on.
    - `cycNumTurnOff`   : Number of cycles of the laser field in the turn-off.
    - `cycNumConst`     : Number of cycles of the laser field in the constant-intensity.
    - `ellip`           : Ellipticity of the laser field [-1≤e≤1, 0 indicates linear polarization and ±1 indicates circular polarization].
    - `azi`             : Azimuth angle of the laser's polarization's principle axis relative to x axis (in radians) (optional, default 0).
    - `t_shift`         : Time shift of the laser (in a.u.) relative to the beginning of TURN-ON (optional, default 0).
    """
    function TrapezoidalLaser(peakInt, waveLen, cycNumTurnOn, cycNumTurnOff, cycNumConst, ellip, azi=0., cep=0., t_shift=0.)
        @assert peakInt>0                                           "[TrapezoidalLaser] Peak intensity must be positive."
        @assert waveLen>0                                           "[TrapezoidalLaser] Wavelength must be positive."
        @assert cycNumTurnOn>0 && cycNumTurnOff>0 && cycNumConst≥0  "[TrapezoidalLaser] Cycle number must be positive."
        @assert -1≤ellip≤1                                          "[TrapezoidalLaser] Ellipticity must be in [-1,1]."
        new(peakInt,waveLen,cycNumTurnOn,cycNumTurnOff,cycNumConst,ellip,azi,cep,t_shift)
    end
    """
    Constructs a new monochromatic elliptically polarized laser field with Trapezoidal-shape envelope.
    # Parameters
    - `peakInt`         : Peak intensity of the laser field (in W/cm²).
    - `WaveLen`         : Wavelength of the laser field (in nm). Must specify either `waveLen` or `angFreq`.
    - `angFreq`         : Angular frequency of the laser field (in a.u.). Must specify either `waveLen` or `angFreq`.
    - `cycNumTurnOn`    : Number of cycles of the laser field in the turn-on.
    - `cycNumTurnOff`   : Number of cycles of the laser field in the turn-off.
    - `cycNumConst`     : Number of cycles of the laser field in the constant-intensity.
    - `ellip`           : Ellipticity of the laser field [-1≤e≤1, 0 indicates linear polarization and ±1 indicates circular polarization].
    - `azi`             : Azimuth angle of the laser's polarization's principle axis relative to x axis (in radians) (optional, default 0).
    - `t_shift`         : Time shift of the laser (in a.u.) relative to the beginning of TURN-ON (optional, default 0).
    """
    function TrapezoidalLaser(; peakInt,
                                waveLen=-1, angFreq=-1,     # must specify either waveLen or angFreq.
                                cycNumTurnOn, cycNumTurnOff, cycNumConst,
                                ellip, azi=0., cep=0., t_shift=0.)
        @assert waveLen>0 || angFreq>0                              "[TrapezoidalLaser] Must specify either waveLen or angFreq."
        @assert cycNumTurnOn>0 && cycNumTurnOff>0 && cycNumConst≥0  "[TrapezoidalLaser] Cycle number must be positive."
        if waveLen>0 && angFreq>0
            @warn "[TrapezoidalLaser] Both waveLen & angFreq are specified, will use waveLen."
        end
        if waveLen==-1
            waveLen = 45.563352525 / angFreq
        end
        TrapezoidalLaser(peakInt,waveLen,cycNumTurnOn,cycNumTurnOff,cycNumConst,ellip,azi,cep,t_shift)
    end
end
"Gets the peak intensity of the laser field (in W/cm²)."
PeakInt(l::TrapezoidalLaser) = l.peak_int
"Gets the wave length of the laser field (in nm)."
WaveLen(l::TrapezoidalLaser) = l.wave_len
"Gets the total cycle number of the laser field."
CycNumTotal(l::TrapezoidalLaser) = l.cyc_num_turn_on + l.cyc_num_turn_off + l.cyc_num_const
"Gets the cycle number of the laser field in the turn-on."
CycNumTurnOn(l::TrapezoidalLaser) = l.cyc_num_turn_on
"Gets the cycle number of the laser field in the turn-off."
CycNumTurnOff(l::TrapezoidalLaser) = l.cycNumTurnOff
"Gets the cycle number of the laser field in the constant-intensity."
CycNumConst(l::TrapezoidalLaser) = l.cycNumConst
"Gets the ellipticity of the laser field."
Ellipticity(l::TrapezoidalLaser) = l.ellip
"Gets the azimuth angle of the laser's polarization's principle axis relative to x axis (in radians)."
Azimuth(l::TrapezoidalLaser) = l.azi
"Gets the angular frequency (ω) of the laser field (in a.u.)."
AngFreq(l::TrapezoidalLaser) = 45.563352525 / l.wave_len
"Gets the period of the laser field (in a.u.)."
Period(l::TrapezoidalLaser) = 2π / AngFreq(l)
"Gets the time shift relative to the beginning of TURN-ON (in a.u.)."
TimeShift(l::TrapezoidalLaser) = l.t_shift
"Gets the peak electric field intensity of the laser field (in a.u.)."
LaserF0(l::TrapezoidalLaser) = sqrt(l.peak_int/(1.0+l.ellip^2)/3.50944521e16)
"Gets the peak vector potential intensity of the laser field (in a.u.)."
LaserA0(l::TrapezoidalLaser) = LaserF0(l) / AngFreq(l)

"Gets the time-dependent x component of the vector potential under dipole approximation."
function LaserAx(l::TrapezoidalLaser)
    local A0 = LaserA0(l); local ω = AngFreq(l); local T = Period(l); local Δt = l.t_shift;
    local N_on = l.cycNumTurnOn; local N_const = l.cycNumConst; local N_off = l.cycNumTurnOff;
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
"Gets the time-dependent y component of the vector potential under dipole approximation."
function LaserAy(l::TrapezoidalLaser)
    local A0 = LaserA0(l); local ω = AngFreq(l); local T = Period(l); local Δt = l.t_shift;
    local N_on = l.cycNumTurnOn; local N_const = l.cycNumConst; local N_off = l.cycNumTurnOff;
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
"Gets the time-dependent x component of the electric field strength under dipole approximation."
function LaserFx(l::TrapezoidalLaser)
    local F0 = LaserF0(l); local ω = AngFreq(l); local T = Period(l); local Δt = l.t_shift;
    local N_on = l.cycNumTurnOn; local N_const = l.cycNumConst; local N_off = l.cycNumTurnOff;
    local t_on = N_on*T; local t_const = N_const*T; local t_off = N_off*T;
    local φ = l.cep; local ε = l.ellip; local ϕ = l.azi;
    return if ϕ==0
        function(t)
            t -= Δt
            F0 * ( ((0<t<t_on)*(t/t_on)+(t_on≤t≤(t_on+t_const))*1.0+((t_on+t_const)<t<(t_on+t_const+t_off))*(1-(t-t_on-t_const)/t_off))*sin(ω*t+φ)*ω - ((0<t<t_on)*(1/t_on)+((t_on+t_const)<t<(t_on+t_const+t_off))*(-1/t_off))*cos(ω*t+φ))
        end
    else
        function(t)
            t -= Δt
            F0 * ( ((0<t<t_on)*(t/t_on)+(t_on≤t≤(t_on+t_const))*1.0+((t_on+t_const)<t<(t_on+t_const+t_off))*(1-(t-t_on-t_const)/t_off))*sin(ω*t+φ)*ω + ((0<t<t_on)*(1/t_on)+((t_on+t_const)<t<(t_on+t_const+t_off))*(-1/t_off))*cos(ω*t+φ)*cos(ϕ) + (((0<t<t_on)*(t/t_on)+(t_on≤t≤(t_on+t_const))*1.0+((t_on+t_const)<t<(t_on+t_const+t_off))*(1-(t-t_on-t_const)/t_off))*cos(ω*t+φ)*ω + ((0<t<t_on)*(1/t_on)+((t_on+t_const)<t<(t_on+t_const+t_off))*(-1/t_off))*sin(ω*t+φ))*-ε*sin(ϕ) )
        end
    end
end
"Gets the time-dependent y component of the electric field strength under dipole approximation."
function LaserFy(l::TrapezoidalLaser)
    local F0 = LaserF0(l); local ω = AngFreq(l); local T = Period(l); local Δt = l.t_shift;
    local N_on = l.cycNumTurnOn; local N_const = l.cycNumConst; local N_off = l.cycNumTurnOff;
    local t_on = N_on*T; local t_const = N_const*T; local t_off = N_off*T;
    local φ = l.cep; local ε = l.ellip; local ϕ = l.azi;
    return if ϕ==0
        function(t)
            t -= Δt
            -F0 * ε * ( ((0<t<t_on)*(t/t_on)+(t_on≤t≤(t_on+t_const))*1.0+((t_on+t_const)<t<(t_on+t_const+t_off))*(1-(t-t_on-t_const)/t_off))*cos(ω*t+φ)*ω + ((0<t<t_on)*(1/t_on)+((t_on+t_const)<t<(t_on+t_const+t_off))*(-1/t_off))*sin(ω*t+φ))
        end
    else
        function(t)
            t -= Δt
            F0 * ( ((0<t<t_on)*(t/t_on)+(t_on≤t≤(t_on+t_const))*1.0+((t_on+t_const)<t<(t_on+t_const+t_off))*(1-(t-t_on-t_const)/t_off))*sin(ω*t+φ)*ω + ((0<t<t_on)*(1/t_on)+((t_on+t_const)<t<(t_on+t_const+t_off))*(-1/t_off))*cos(ω*t+φ)*-sin(ϕ) + (((0<t<t_on)*(t/t_on)+(t_on≤t≤(t_on+t_const))*1.0+((t_on+t_const)<t<(t_on+t_const+t_off))*(1-(t-t_on-t_const)/t_off))*cos(ω*t+φ)*ω + ((0<t<t_on)*(1/t_on)+((t_on+t_const)<t<(t_on+t_const+t_off))*(-1/t_off))*sin(ω*t+φ))*-ε*cos(ϕ) )
        end
    end
end

"Prints the information about the laser."
Base.show(io::IO, l::TrapezoidalLaser) = print(io,"[MonochromaticLaser] Envelope Trapezoidal, Wavelength=$(l.waveLen) nm, TurnOn/Constant/TurnOff: $(l.cycNumTurnOn)/$(l.cycNumConst)/$(l.cycNumTurnOff) cycle(s), e=$(l.ellip)"
                                                * (l.ellip==0 ? " [Linearly polarized]" : "") * (abs(l.ellip)==1 ? " [Circularly polarized]" : "")
                                                * ", PrincipleAxisAzimuth=$(l.azi/π*180)°" * (l.t_shift==0 ? "" : ", Rises at t₀=$(l.t_shift) a.u.") * (l.cep==0 ? "" : ", CEP=$(l.cep)"))
