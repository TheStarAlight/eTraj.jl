
"""
    module Lasers

The `Laser` module contains abstraction of laser and provides some pre-defined lasers for use.
"""
module Lasers

include("imports.jl")
include("exports.jl")

"""
    abstract type Laser

Represents an abstract laser, supertype of all lasers.
"""
abstract type Laser end
"""
    abstract type MonochromaticLaser <: Laser

Represents an abstract monochromatic laser.
"""
abstract type MonochromaticLaser <: Laser end

include("Cos4Laser.jl")
include("Cos2Laser.jl")
include("TrapezoidalLaser.jl")
include("GaussianLaser.jl")
include("BichromaticLaser.jl")
include("docs.jl")

end