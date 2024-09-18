
# ==== General Laser ====

@doc """
    LaserAx(l::Laser) -> Ax(t)

Gets the time-dependent x component of the vector potential under dipole approximation.
"""
LaserAx

@doc """
    LaserAy(l::Laser) -> Ay(t)

Gets the time-dependent y component of the vector potential under dipole approximation.
"""
LaserAy

@doc """
    LaserFx(l::Laser) -> Fx(t)

Gets the time-dependent x component of the electric field strength under dipole approximation.
"""
LaserFx

@doc """
    LaserFy(l::Laser) -> Fy(t)

Gets the time-dependent y component of the electric field strength under dipole approximation.
"""
LaserFy

# ==== MonochromaticLaser ====

@doc """
    PeakInt(l::MonochromaticLaser)

Gets the peak intensity of the laser (in W/cm²).
"""
PeakInt

@doc """
    WaveLen(l::MonochromaticLaser)

Gets the wave length of the laser (in nm).
"""
WaveLen

@doc """
    LaserF0(l::MonochromaticLaser)

Gets the peak electric field intensity of the laser (in a.u.).
"""
LaserF0

@doc """
    LaserA0(l::MonochromaticLaser)

Gets the peak vector potential intensity of the laser (in a.u.).
"""
LaserA0

@doc """
    Ellipticity(l::MonochromaticLaser)

Gets the ellipticity of the laser.
"""
Ellipticity

@doc """
    AngFreq(l::MonochromaticLaser)

Gets the angular frequency (ω) of the laser (in a.u.).
"""
AngFreq

@doc """
    Period(l::MonochromaticLaser)

Gets the period of the laser (in a.u.).
"""
Period

@doc """
    CEP(l::MonochromaticLaser)

Gets the Carrier-Envelope Phase (CEP) of the laser.
"""
CEP

@doc """
    Azimuth(l::MonochromaticLaser)

Gets the azimuth angle of the laser's polarization's principle axis relative to x axis (in radians).
"""
Azimuth

@doc """
    KeldyshParameter(l::MonochromaticLaser, Ip)

Gets the Keldysh parameter γ₀ of the laser, given the ionization energy `Ip` (in a.u.).
"""
KeldyshParameter

@doc """
    UnitEnvelope(l::MonochromaticLaser) -> env(t)

Gets the unit envelope function (the peak value is 1) of the laser field.
"""
UnitEnvelope

# ==== Cos2, Cos4 and Gaussian ====

@doc """
    TimeShift(l::{Cos2Laser, Cos4Laser, GaussianLaser})

Gets the time shift of the laser relative to the peak (in a.u.).
"""
TimeShift

# ==== Cos2, Cos4 ====

@doc """
    CycNum(l::{Cos2Laser, Cos4Laser})

Gets the cycle number of the laser.
"""
CycNum

# ==== Trapezoidal ====

@doc """
    CycNumTurnOn(l::TrapezoidalLaser)

Gets the cycle number of the laser in the turn-on.
"""
CycNumTurnOn

@doc """
    CycNumConst(l::TrapezoidalLaser)

Gets the cycle number of the laser in the constant-intensity stage.
"""
CycNumConst

@doc """
    CycNumTurnOff(l::TrapezoidalLaser)

Gets the cycle number of the laser in the turn-off.
"""
CycNumTurnOff

@doc """
    CycNumTotal(l::TrapezoidalLaser)

Gets the total cycle number of the laser.
"""
CycNumTotal

@doc """
    TimeTurnOn(l::TrapezoidalLaser)

Gets the time shift relative to the beginning of turn-on (in a.u.).
"""
TimeTurnOn

# ==== Gaussian ====

@doc """
    SpreadCycNum(l::GaussianLaser)

Gets the half-temporal width of the laser, namely σ (converted to cycle numbers).
"""
SpreadCycNum

@doc """
    SpreadDuration(l::GaussianLaser)

Gets the half-temporal width of the laser, namely σ (in a.u.).
"""
SpreadDuration

@doc """
    FWHM_Duration(l::GaussianLaser)

Gets the temporal FWHM (Full Width at Half Maxima) of the laser (in a.u.).
"""
FWHM_Duration

# ==== BichromaticLaser ====

@doc """
    Laser1(l::BichromaticLaser)

Gets the first `MonochromaticLaser` component in `l`.
"""
Laser1

@doc """
    Laser2(l::BichromaticLaser)

Gets the second `MonochromaticLaser` component in `l`.
"""
Laser2

@doc """
    Delay21(l::BichromaticLaser)

Gets the delay of the second laser respective to the first in `l`.
"""
Delay21
