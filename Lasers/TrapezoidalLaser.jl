
"Represents a monochromatic elliptically polarized laser field with Trapezoidal-shape envelope propagating in z direction."
struct TrapezoidalLaser <: MonochromaticLaser
    "Peak intensity of the laser field (in W/cm^2)."
    peakInt;
    "Wavelength of the laser field (in NANOMETER)."
    waveLen;
    "Cycle number of the laser field in the turn-on."
    cycNumTurnOn;
    "Cycle number of the laser field in the turn-off."
    cycNumTurnOff;
    "Cycle number of the laser field in the constant-intensity."
    cycNumConst;
    "Ellpticity of the laser field."
    ellip;
    "Carrier-Envelope-Phase (CEP) of the laser field."
    cep;
    """
    Constructs a new monochromatic elliptically polarized laser field with Trapezoidal-shape envelope.
    # Parameters
    - `peakInt`         : Peak intensity of the laser field (in W/cm²).
    - `WaveLen`         : Wavelength of the laser field (in nm).
    - `cycNumTurnOn`    : Number of cycles of the laser field in the turn-on.
    - `cycNumTurnOff`   : Number of cycles of the laser field in the turn-off.
    - `cycNumConst`     : Number of cycles of the laser field in the constant-intensity.
    - `Ellpticity`      : Ellpticity of the laser field [0≤e≤1, 0 indicates linear polarization (in x direction) and 1 indicates circular polarization].
    - `cep`             : Carrier-Envelope-Phase of the laser field (optional, default 0).
    """
    function TrapezoidalLaser(peakInt, waveLen, cycNumTurnOn, cycNumTurnOff, cycNumConst, ellip, cep=0.)
        @assert peakInt>0   "[TrapezoidalLaser] Peak intensity must be positive."
        @assert waveLen>0   "[TrapezoidalLaser] Wavelength must be positive."
        @assert cycNumTurnOn>0 && cycNumTurnOff>0 && cycNumConst≥0  "[TrapezoidalLaser] Cycle number must be positive."
        @assert 0≤ellip≤1   "[TrapezoidalLaser] Ellpticity must be in [0,1]."
        new(peakInt,waveLen,cycNumTurnOn,cycNumTurnOff,cycNumConst,ellip,cep)
    end
    function TrapezoidalLaser(;peakInt, waveLen, cycNumTurnOn, cycNumTurnOff, cycNumConst, ellip, cep=0.)
        TrapezoidalLaser(peakInt,waveLen,cycNumTurnOn,cycNumTurnOff,cycNumConst,ellip,cep)
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
    - `Ellpticity`      : Ellpticity of the laser field [0≤e≤1, 0 indicates linear polarization (in x direction) and 1 indicates circular polarization].
    - `cep`             : Carrier-Envelope-Phase of the laser field (optional, default 0).
    """
    function TrapezoidalLaser(; peakInt,
                                waveLen=-1, angFreq=-1,     # must specify either waveLen or angFreq.
                                cycNumTurnOn, cycNumTurnOff, cycNumConst, ellip,
                                cep=0.)
        @assert peakInt>0               "[TrapezoidalLaser] Peak intensity must be positive."
        @assert waveLen>0 || angFreq>0  "[TrapezoidalLaser] Must specify either waveLen or angFreq."
        @assert cycNumTurnOn>0 && cycNumTurnOff>0 && cycNumConst≥0  "[TrapezoidalLaser] Cycle number must be positive."
        @assert 0≤ellip≤1               "[TrapezoidalLaser] Ellpticity must be in [0,1]."
        if waveLen>0 && angFreq>0
            @warn "[TrapezoidalLaser] Both waveLen & angFreq are specified, will use waveLen."
        end
        if waveLen==-1
            waveLen = 45.563352525 / angFreq
        end
        new(peakInt,waveLen,cycNumTurnOn,cycNumTurnOff,cycNumConst,ellip,cep)
    end
end
"Gets the peak intensity of the laser field (in W/cm²)."
PeakInt(l::TrapezoidalLaser) = l.peakInt
"Gets the wave length of the laser field (in nm)."
WaveLen(l::TrapezoidalLaser) = l.waveLen
"Gets the total cycle number of the laser field."
CycNumTotal(l::TrapezoidalLaser) = l.cycNumTurnOn + l.cycNumTurnOff + l.cycNumConst
"Gets the cycle number of the laser field in the turn-on."
CycNumTurnOn(l::TrapezoidalLaser) = l.cycNumTurnOn
"Gets the cycle number of the laser field in the turn-off."
CycNumTurnOff(l::TrapezoidalLaser) = l.cycNumTurnOff
"Gets the cycle number of the laser field in the constant-intensity."
CycNumConst(l::TrapezoidalLaser) = l.cycNumConst
"Gets the ellpticity of the laser field."
Ellpticity(l::TrapezoidalLaser) = l.ellip
"Gets the angular frequency (ω) of the laser field (in a.u.)."
AngFreq(l::TrapezoidalLaser) = 45.563352525 / l.waveLen
"Gets the period of the laser field (in a.u.)."
Period(l::TrapezoidalLaser) = 2π / AngFreq(l)
"Gets the peak electric field intensity of the laser field (in a.u.)."
LaserF0(l::TrapezoidalLaser) = sqrt(l.peakInt/(1.0+l.ellip^2)/3.50944521e16)
"Gets the peak vector potential intensity of the laser field (in a.u.)."
LaserA0(l::TrapezoidalLaser) = LaserF0(l) / AngFreq(l)

"Gets the time-dependent x component of the vector potential under dipole approximation."
function LaserAx(l::TrapezoidalLaser)
    local A0 = LaserA0(l); local ω = AngFreq(l); local T = Period(l);
    local N_on = l.cycNumTurnOn; local N_const = l.cycNumConst; local N_off = l.cycNumTurnOff;
    local t_on = N_on*T; local t_const = N_const*T; local t_off = N_off*T;
    local φ = l.cep;
    return function(t)  # this function might be accessed by GPU, thus the conditional exp. should be avoided.
        A0 * ( (0<t<t_on)*(t/t_on) + (t_on≤t≤(t_on+t_const))*1.0 + ((t_on+t_const)<t<(t_on+t_const+t_off))*(1-(t-t_on-t_const)/t_off) ) * cos(ω*t+φ)
    end
end
"Gets the time-dependent y component of the vector potential under dipole approximation."
function LaserAy(l::TrapezoidalLaser)
    local A0 = LaserA0(l); local ω = AngFreq(l); local T = Period(l);
    local N_on = l.cycNumTurnOn; local N_const = l.cycNumConst; local N_off = l.cycNumTurnOff;
    local t_on = N_on*T; local t_const = N_const*T; local t_off = N_off*T;
    local φ = l.cep; local ε = l.ellip;
    return function(t)
        A0 * ( (0<t<t_on)*(t/t_on) + (t_on≤t≤(t_on+t_const))*1.0 + ((t_on+t_const)<t<(t_on+t_const+t_off))*(1-(t-t_on-t_const)/t_off) ) * sin(ω*t+φ) * ε
    end
end
"Gets the time-dependent x component of the electric field strength under dipole approximation."
function LaserFx(l::TrapezoidalLaser)
    local F0 = LaserF0(l); local ω = AngFreq(l); local T = Period(l);
    local N_on = l.cycNumTurnOn; local N_const = l.cycNumConst; local N_off = l.cycNumTurnOff;
    local t_on = N_on*T; local t_const = N_const*T; local t_off = N_off*T;
    local φ = l.cep;
    return function(t)
        F0 * ( ((0<t<t_on)*(t/t_on)+(t_on≤t≤(t_on+t_const))*1.0+((t_on+t_const)<t<(t_on+t_const+t_off))*(1-(t-t_on-t_const)/t_off))*sin(ω*t+φ)*ω - ((0<t<t_on)*(1/t_on)+((t_on+t_const)<t<(t_on+t_const+t_off))*(-1/t_off))*cos(ω*t+φ))
    end
end
"Gets the time-dependent y component of the electric field strength under dipole approximation."
function LaserFy(l::TrapezoidalLaser)
    local F0 = LaserF0(l); local ω = AngFreq(l); local T = Period(l);
    local N_on = l.cycNumTurnOn; local N_const = l.cycNumConst; local N_off = l.cycNumTurnOff;
    local t_on = N_on*T; local t_const = N_const*T; local t_off = N_off*T;
    local φ = l.cep; local ε = l.ellip;
    return function(t)
        -F0 * ε * ( ((0<t<t_on)*(t/t_on)+(t_on≤t≤(t_on+t_const))*1.0+((t_on+t_const)<t<(t_on+t_const+t_off))*(1-(t-t_on-t_const)/t_off))*cos(ω*t+φ)*ω + ((0<t<t_on)*(1/t_on)+((t_on+t_const)<t<(t_on+t_const+t_off))*(-1/t_off))*sin(ω*t+φ))
    end
end

"Prints the information about the laser."
Base.show(io::IO, l::TrapezoidalLaser) = print(io,"[MonochromaticLaser] Envelope Trapezoidal, Wavelength=$(l.waveLen) nm, TurnOn/Constant/TurnOff: $(l.cycNumTurnOn)/$(l.cycNumConst)/$(l.cycNumTurnOff) cycle(s), e=$(l.ellip)"
                                               * (l.ellip==0 ? " [Linearly (x ax.) polarized]" : "") * (l.ellip==1 ? " [Circularly polarized]" : "") * (l.cep==0 ? "" : ", CEP=$(l.cep)"))