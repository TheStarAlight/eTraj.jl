

"""
The Targets module provides information about the atoms or molecules.
"""
module Targets

include("imports.jl")
include("exports.jl")

abstract type Target end
abstract type SAEAtomBase <: Target end
abstract type MoleculeBase <: Target end

include("SAEAtom.jl")
include("HydrogenLikeAtom.jl")
include("SAEAtomBase_shared.jl")
include("AtomDatabase.jl")

include("MolecularCalculators/MolecularCalculatorBase.jl")
include("GenericMolecule.jl")

end