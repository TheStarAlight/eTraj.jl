

"""
The Targets module provides information about the targeting atoms or molecules.
"""
module Targets

export Target, SAEAtom, HydrogenLikeAtom
export IonPotential, AsympNuclCharge, TargetPotential, TargetForce, ADKRateExp


abstract type Target end


"Represents an atom under single-active-electron (SAE) approximation."
abstract type SAEAtom <: Target end
# should implement TargetPotential, TargetForce, ADKRateExp.

include("HydrogenLikeAtom.jl")

end