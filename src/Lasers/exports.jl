
# Laser types
export Laser, MonochromaticLaser, BichromaticLaser
export Cos4Laser, Cos2Laser, GaussianLaser
# General properties
export LaserAx, LaserAy, LaserFx, LaserFy
export Serialize
# Laser-specific properties
# MonochromaticLaser
export LaserF0, LaserA0, PeakInt, WaveLen, Ellipticity, AngFreq, Period, CEP, Azimuth, KeldyshParameter, UnitEnvelope
# Cos2, Cos4, Gaussian
export TimeShift
# Cos2, Cos4
export CycNum
# Gaussian
export SpreadCycNum, SpreadDuration, FWHM_Duration
# BichromaticLaser
export Laser1, Laser2, Delay21