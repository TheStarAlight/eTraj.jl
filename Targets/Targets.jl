

"""
The Targets module provides information about the targeting atoms or molecules.
"""
module Targets

export Target, SAEAtom, HydrogenLikeAtom
export IonPotential, AsympNuclCharge, TargetPotential, TargetForce, TrajectoryFunction, ADKRateExp


abstract type Target end


"Represents an atom under single-active-electron (SAE) approximation."
abstract type SAEAtom <: Target end
# should implement TargetPotential, TargetForce, TrajectoryFunction, ADKRateExp.

include("HydrogenLikeAtom.jl")

end