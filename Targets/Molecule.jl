
"Represents a molecule."
struct Molecule <: Target
    "Atoms in the molecule with their coordinates, which are stored as a vector of (:AtomName, X, Y, Z)."
    atoms;
    "Total charge of the molecule (ion)."
    molCharge;
    "Euler angles (ZYZ convention) specifying the molecule's orientation."
    rot_α;
    rot_β;
    rot_γ;
    "Name of the molecule."
    name::String;

    """
    Initializes a new instance of `Molecule` with given parameters.
    # Parameters
    - `atoms`                   : Atoms in the molecule with their coordinates, which are stored as a vector of (:AtomName, X, Y, Z).
    - `molCharge`               : Total charge of the molecule (ion) (optional, default 0).
    - `rot_α`,`rot_β`,`rot_γ`   : Euler angles (ZYZ convention) specifying the molecule's orientation (optional, default 0).
    - `name`                    : Name of the molecule (optional).
    """
    function Molecule(atoms,molCharge=0,rot_α=0.,rot_β=0.,rot_γ=0.,name::String="[NA]")
        new(atoms,molCharge,rot_α,rot_β,rot_γ,name)
    end
    function Molecule(;atoms,molCharge=0,rot_α=0.,rot_β=0.,rot_γ=0.,name::String="[NA]")
        new(atoms,molCharge,rot_α,rot_β,rot_γ,name)
    end
    function Molecule(;atoms,molCharge=0,rot::Tuple=(0.,0.,0.),name::String="[NA]")
        new(atoms,molCharge,rot[1],rot[2],rot[3],name)
    end
end

"Gets the atoms in the molecule with their coordinates."
Atoms(mol::Molecule) = mol.atoms
"Gets the total charge of the molecule (ion)."
MolCharge(mol::Molecule) = mol.molCharge
"Gets the Euler angles (ZYZ convention) specifying the molecule's orientation in format (α,β,γ)."
Rotation(mol::Molecule) = (mol.rot_α,mol.rot_β,mol.rot_γ)
"Gets the name of the molecule."
TargetName(mol::Molecule) = mol.name

"""
Exports the given molecule's atom information to string as `MolecularCalculator`'s input.
Note: Rotations defined by the Euler angles wouldn't be applied.
"""
function exportMolAtomInfo(mol::Molecule)
    atomToString(atom) = join([String(atom[1]),atom[2],atom[3],atom[4]], " ")
    return join(map(atomToString, mol.atoms),";")
end

#TODO: The structure of Molecule needs to be improved.
#TODO: Maybe init the Molecule with a precomputed data is a better choice.
#TODO: If there's no data, the user may provide information about this molecule, and the MolecularCalculator would be utilized to calculate.
#TODO: And what to calculate can be specified. The user can append data by loading the Molecule data and call methods in the MolecularCalculator.
