

"""
The Targets module provides information about the targeting atoms or molecules.
"""
module Targets

export Target, SAEAtomBase
export SAEAtom, HydrogenLikeAtom
export IonPotential, AsympNuclCharge, TargetPotential, TargetForce, TrajectoryFunction, ADKRateExp

abstract type Target end
abstract type SAEAtomBase <: Target end

include("SAEAtom.jl")
include("HydrogenLikeAtom.jl")

end