
"""
The Laser module provides information about laser fields.
"""
module Lasers

export Laser, Cos4Laser
export PeakInt, WaveLen, CycNum, Ellpticity, AngFreq, Period
export LaserF0, LaserA0, LaserAx, LaserAy, LaserFx, LaserFy


abstract type Laser end
abstract type MonochromaticLaser <: Laser end
abstract type BichromaticLaser <: Laser end

include("Cos4Laser.jl")

end