
"""
    get_mol(name::String; [rot_α=0.0] [,rot_β=0.0] [,rot_γ=0.0])

Constructs a `GenericMolecule` in the database based on the `name` given.
Euler angles which describe the molecule's rotation can be applied via the keyword arguments `rot_α`, `rot_β` and `rot_γ`.
The available molecule names can be found by invoking [`get_available_mols`](@ref).

## Examples

```jldoctest

julia> get_mol("Hydrogen")
[GenericMolecule] Hydrogen (H₂)
Asymp coeff of HOMO available
WFAT data of HOMO available
#          E (Ha)  occp
⋮  ⋮         ⋮      ⋮⋮
3  LUMO+1   0.232  ----
2  LUMO     0.148  ----
1  HOMO    -0.594  -↿⇂-

julia> get_mol("Oxygen"; rot_β=π)
[GenericMolecule] Oxygen (O₂), αβγ=(0.0°,180.0°,0.0°)
Asymp coeff of α-HOMO-1 & α-HOMO available
WFAT data of α-HOMO-1 & α-HOMO available
#          Eα(Ha)  occp  Eβ(Ha)
⋮    ⋮        ⋮     ⋮⋮      ⋮     ⋮
11 LUMO+1   0.681  ----   0.741 LUMO+3
10 LUMO     0.377  ----   0.431 LUMO+2
9  HOMO    -0.554  -↿--   0.110 LUMO+1
8  HOMO-1  -0.554  -↿--   0.110 LUMO
7  HOMO-2  -0.763  -↿⇂-  -0.575 HOMO
6  HOMO-3  -0.841  -↿⇂-  -0.575 HOMO-1
5  HOMO-4  -0.841  -↿⇂-  -0.702 HOMO-2
4  HOMO-5  -1.204  -↿⇂-  -0.993 HOMO-3
⋮    ⋮        ⋮     ⋮⋮      ⋮     ⋮
```

See also [`get_available_mols`](@ref) and [`GenericMolecule`](@ref).
"""
function get_mol(name::String; rot_α=0.,rot_β=0.,rot_γ=0.)
    if name in keys(MoleculeDataPath)
        return LoadMolecule(join([Targets|>Base.moduleroot|>pathof|>dirname, "/Targets/Database/", MoleculeDataPath[name]]); rot_α,rot_β,rot_γ)
    end
    error("[get_mol] Molecule name `$name` not found in database.")
end

function get_available_mols()
    return collect(keys(MoleculeDataPath))
end

const MoleculeDataPath = OrderedDict{String, String}(
    "Hydrogen"          => "Molecule_Hydrogen.jld2",
    "Nitrogen"          => "Molecule_Nitrogen.jld2",
    "Oxygen"            => "Molecule_Oxygen.jld2",
    "Carbon Monoxide"   => "Molecule_CarbonMonoxide.jld2",
    "Nitric Oxide"      => "Molecule_NitricOxide.jld2",
    "Hydrochloric Acid" => "Molecule_HydrochloricAcid.jld2",
    "Carbon Dioxide"    => "Molecule_CarbonDioxide.jld2",
    "Sulfur Dioxide"    => "Molecule_SulfurDioxide.jld2",
    "Water"             => "Molecule_Water.jld2",
    "Ammonia"           => "Molecule_Ammonia.jld2",
    "Acetylene"         => "Molecule_Acetylene.jld2",
    "Methane"           => "Molecule_Methane.jld2",
    "Benzene"           => "Molecule_Benzene.jld2"
)
