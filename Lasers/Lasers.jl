
"""
The Laser module provides information about laser fields.
"""
module Lasers

export Laser, Cos4Laser
export PeakInt, WaveLenNM, CycNum, Ellpticity, AngFreq, Period
export LaserF0, LaserA0, LaserAx, LaserAy, LaserFx, LaserFy


"Represents an abstract laser."
abstract type Laser end

include("Cos4Laser.jl")

end