
"Represents a monochromatic elliptically polarized laser field with Cos4-shape envelope propagating in z direction."
struct Cos4Laser <: MonochromaticLaser
    "Peak intensity of the laser field (in W/cm^2)."
    peakInt;
    "Wavelength of the laser field (in NANOMETER)."
    waveLen;
    "Cycle number of the laser field."
    cycNum;
    "Ellpticity of the laser field."
    ellip;
    "Carrier-Envelope-Phase (CEP) of the laser field."
    cep;
    "Constructs a new monochromatic elliptically polarized laser field."
    Cos4Laser(peakInt, waveLen, cycNum, ellip, cep) = new(peakInt,waveLen,cycNum,ellip,cep);
    #TODO: add support for more flexible constructor function.
end
"Gets the peak intensity of the laser field (in W/cm^2)."
PeakInt(l::Cos4Laser) = l.peakInt
"Gets the wave length of the laser field (in nm)."
WaveLen(l::Cos4Laser) = l.waveLen
"Gets the cycle number of the laser field."
CycNum(l::Cos4Laser) = l.cycNum
"Gets the ellpticity of the laser field."
Ellpticity(l::Cos4Laser) = l.ellip
"Gets the angular frequency (ω) of the laser field (in a.u.)."
AngFreq(l::Cos4Laser) = 45.563352525 / l.waveLen
"Gets the period of the laser field (in a.u.)."
Period(l::Cos4Laser) = 2π / AngFreq(l)
"Gets the peak electric field intensity of the laser field (in a.u.)."
LaserF0(l::Cos4Laser) = sqrt(l.peakInt/(1.0+l.ellip^2)/3.50944521e16)
"Gets the peak vector potential intensity of the laser field (in a.u.)."
LaserA0(l::Cos4Laser) = LaserF0(l) / AngFreq(l)
#TODO: make further consideration on providing necessary property accessors.

"Gets the time-dependent x component of the vector potential under dipole approximation."
function LaserAx(l::Cos4Laser)
    local A0 = LaserA0(l); local ω = AngFreq(l); local N = l.cycNum; local φ = l.cep;
    return function(t)
        A0 * cos(ω*t/(2N))^4 * (abs(ω*real(t))<N*π) * cos(ω*t+φ)
    end
end
"Gets the time-dependent y component of the vector potential under dipole approximation."
function LaserAy(l::Cos4Laser)
    local A0 = LaserA0(l); local ω = AngFreq(l); local N = l.cycNum; local φ = l.cep; local ε = l.ellip;
    return function(t)
        A0 * cos(ω*t/(2N))^4 * (abs(ω*real(t))<N*π) * sin(ω*t+φ) * ε
    end
end
"Gets the time-dependent x component of the electric field strength under dipole approximation."
function LaserFx(l::Cos4Laser)
    local F0 = LaserF0(l); local ω = AngFreq(l); local N = l.cycNum; local φ = l.cep;
    return function(t)
        F0 * cos(ω*t/(2N))^3 * (abs(ω*real(t))<N*π) * ( cos(ω*t/(2N))*sin(ω*t+φ) + 2/N*sin(ω*t/(2N))*cos(ω*t+φ))
    end
end
"Gets the time-dependent y component of the electric field strength under dipole approximation."
function LaserFy(l::Cos4Laser)
    local F0 = LaserF0(l); local ω = AngFreq(l); local N = l.cycNum; local φ = l.cep; local ε = l.ellip;
    return function(t)
        F0 * cos(ω*t/(2N))^3 * (abs(ω*real(t))<N*π) * (-cos(ω*t/(2N))*cos(ω*t+φ) + 2/N*sin(ω*t/(2N))*sin(ω*t+φ)) * ε
    end
end