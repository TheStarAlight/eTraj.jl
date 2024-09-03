
# Laser types
export Laser, MonochromaticLaser, BichromaticLaser
export Cos4Laser, Cos2Laser, TrapezoidalLaser, GaussianLaser
# General properties
export LaserAx, LaserAy, LaserFx, LaserFy
export Serialize
# Laser-specific properties
# MonochromaticLaser
export LaserF0, LaserA0, PeakInt, WaveLen, Ellipticity, AngFreq, Period, CEP, Azimuth, TimeShift, KeldyshParameter
export UnitEnvelope
# Cos2, Cos4
export CycNum
# Trapezoidal
export CycNumTurnOn, CycNumConst, CycNumTurnOff, CycNumTotal
# Gaussian
export SpreadCycNum, SpreadDuration, FWHM_Duration
# BichromaticLaser
export Laser1, Laser2, Delay21