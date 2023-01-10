
"Represents a molecule."
struct Molecule <: Target
    "Atoms in the molecule, which are stored as a vector of (:AtomName, X, Y, Z)."
    atoms;
    name::String;
end