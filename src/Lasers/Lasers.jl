
"""
The Laser module provides information about laser fields.
"""
module Lasers

include("imports.jl")
include("exports.jl")

abstract type Laser end
abstract type MonochromaticLaser <: Laser end

include("Cos4Laser.jl")
include("Cos2Laser.jl")
include("TrapezoidalLaser.jl")
include("GaussianLaser.jl")
include("BichromaticLaser.jl")

end