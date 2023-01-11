

"""
The Targets module provides information about the targeting atoms or molecules.
"""
module Targets

export Target, SAEAtomBase
export SAEAtom, HydrogenLikeAtom
export HAtom, He1pAtom, Li2pAtom, HeAtom, NeAtom, Ne1pAtom, Ne2pAtom, ArAtom, Ar1pAtom, Ar2pAtom, VAtom, NiAtom, KrAtom, Kr1pAtom, RbAtom, NbAtom, PdAtom, XeAtom, Xe1pAtom, TaAtom
export IonPotential, AsympNuclCharge, TargetName, TargetPotential, TargetForce, TrajectoryFunction, ADKRateExp

export Molecule
export Atoms, MolCharge, exportMolAtomInfo

abstract type Target end
abstract type SAEAtomBase <: Target end

include("SAEAtom.jl")
include("HydrogenLikeAtom.jl")
include("AtomLibrary.jl")
include("Molecule.jl")

end