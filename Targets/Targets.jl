

"""
The Targets module provides information about the targeting atoms or molecules.
"""
module Targets

export Target
export IonPotential, AsympNuclCharge, TargetName, TargetPotential, TargetForce, TrajectoryFunction

export SAEAtomBase, SAEAtom, HydrogenLikeAtom
export HAtom, He1pAtom, Li2pAtom, HeAtom, NeAtom, Ne1pAtom, Ne2pAtom, ArAtom, Ar1pAtom, Ar2pAtom, VAtom, NiAtom, KrAtom, Kr1pAtom, RbAtom, NbAtom, PdAtom, XeAtom, Xe1pAtom, TaAtom
export ADKRateExp

export Molecule
export MolAtoms, MolAtomCoords, MolCharge, MolEnergyLevels, MolHOMOEnergy, MolHOMOIndex, MolWFATAvailableIndices, MolWFATData, MolWFATStructureFactor_G, MolRotation, SetMolRotation, MolExportAtomInfo
export MolCalcEnergyData!, MolCalcWFATData!, MolSaveDataAs

export MolecularCalculators
export PySCFMolecularCalculator

abstract type Target end
abstract type SAEAtomBase <: Target end

include("SAEAtom.jl")
include("HydrogenLikeAtom.jl")
include("AtomLibrary.jl")

include("MolecularCalculators/MolecularCalculator.jl")
include("Molecule.jl")

end