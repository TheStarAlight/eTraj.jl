

"""
    module Targets

The `Targets` module contains abstraction of the targets and provides some pre-defined atoms or molecules.
"""
module Targets

include("imports.jl")
include("exports.jl")

"""
    abstract type Target

Represents an abstract target, supertype of all targets.
"""
abstract type Target end
"""
    abstract type SAEAtomBase <: Target

Represents an abstract atom under single-active-electron approximation.
"""
abstract type SAEAtomBase <: Target end
"""
    abstract type MoleculeBase <: Target

Represents an abstract molecule.
"""
abstract type MoleculeBase <: Target end

include("SAEAtom.jl")
include("HydrogenLikeAtom.jl")
include("SAEAtomBase_shared.jl")
include("Database/AtomDatabase.jl")
include("Database/MoleculeDatabase.jl")

include("MolecularCalculators/MolecularCalculatorBase.jl")
include("GenericMolecule.jl")

include("docs.jl")

end