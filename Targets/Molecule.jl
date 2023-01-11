
"Represents a molecule."
struct Molecule <: Target
    "Atoms in the molecule with their coordinates, which are stored as a vector of (:AtomName, X, Y, Z)."
    atoms;
    "Total charge of the molecule (ion)."
    molCharge;
    "Name of the molecule."
    name::String;

    """
    Initializes a new instance of `Molecule` with given parameters.
    # Parameters
    - `atoms`       : Atoms in the molecule with their coordinates, which are stored as a vector of (:AtomName, X, Y, Z).
    - `molCharge`   : Total charge of the molecule (ion) (optional, default 0).
    - `name`        : Name of the molecule (optional).
    """
    function Molecule(atoms,molCharge=0,name::String="[NA]")
        new(atoms,molCharge,name)
    end
    function Molecule(;atoms,molCharge=0,name::String="[NA]")
        new(atoms,molCharge,name)
    end
end

"Gets the atoms in the molecule with their coordinates."
Atoms(mol::Molecule) = mol.atoms
"Gets the total charge of the molecule (ion)."
MolCharge(mol::Molecule) = mol.molCharge
"Gets the name of the molecule."
TargetName(mol::Molecule) = mol.name

"Exports the given molecule's atom information to string as `MolecularCalculator`'s input."
function exportMolAtomInfo(mol::Molecule)
    atomToString(atom) = join([String(atom[1]),atom[2],atom[3],atom[4]], " ")
    return join(map(atomToString, mol.atoms),";")
end