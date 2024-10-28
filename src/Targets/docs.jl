
# ==== General Target ====

@doc """
    IonPotential(t::SAEAtomBase)

Gets the ionization potential of the atom.

    IonPotential(mol::GenericMolecule)

Gets the ionization potential of the molecule's HOMO.

    IonPotential(mol::GenericMolecule, orbit_ridx)

Gets the ionization potential of the specified MO of molecule.
- `orbit_ridx`: Index of selected orbit relative to the HOMO (e.g., `0` indicates HOMO, and `-1` indicates HOMO-1).
                For open-shell molecules, according to α/β spins, should be passed in format `(spin, idx)` where for α orbitals `spin=1` and for β orbitals `spin=2`.
"""
IonPotential

@doc """
    AsympNuclCharge(t::Target)

Gets the asymptotic nuclear charge of the target (after an electron got ionized).
"""
AsympNuclCharge

@doc """
    TargetName(t::Target)

Gets the name of the target.
"""
TargetName

@doc """
    TargetPotential(t::SAEAtomBase) -> V(x,y,z)

Gets the potential function of the atom.

    TargetPotential(mol::GenericMolecule) -> V(x,y,z)

Gets the asymptotic Coulomb potential function of the molecule.
"""
TargetPotential

@doc """
    TargetForce(t::SAEAtomBase) -> F(x,y,z) -> (Fx,Fy,Fz)

Gets the Coulomb force exerted on the electron from the atom.

    TargetForce(mol::GenericMolecule) -> F(x,y,z) -> (Fx,Fy,Fz)

Gets the asymptotic Coulomb force exerted on the electron from the molecular ion.
"""
TargetForce

@doc """
    TrajectoryFunction(t::Target, dimension = 2|3, laserFx::Function, laserFy::Function, phase_method = :CTMC|:QTMC|:SCTS)

Gets the trajectory function of the electron according to given parameter.
"""
TrajectoryFunction

# ==== SAEAtomBase ====

@doc """
    AngularQuantumNumber(t::SAEAtomBase)

Gets the angular quantum number (l) of the atom.
"""
AngularQuantumNumber

@doc """
    MagneticQuantumNumber(t::SAEAtomBase)

Gets the magnetic quantum number (m) of the atom.
"""
MagneticQuantumNumber

@doc """
    AsympCoeff(t::SAEAtomBase)

Gets the asymptotic coefficient (C_κl) of the atom.
"""
AsympCoeff

@doc """
    SoftCore(t::SAEAtomBase)

Gets the soft core parameter of the atom.
"""
SoftCore

@doc """
    QuantizationAxisOrientaion(t::SAEAtomBase)

Gets the orientation of the quantization axis of the atom in spherical coordinates (θ,ϕ).
"""
QuantizationAxisOrientaion

# ==== GenericMolecule ====

@doc """
    MolAtoms(mol::GenericMolecule) -> ["atom1","atom2",...]

Gets the atoms in the molecule with their coordinates.
"""
MolAtoms

@doc """
    MolAtomCoords(mol::GenericMolecule) -> [x1 y1 z1; x2 y2 z2; ...]

Gets the atoms' coordinates in the molecule.
"""
MolAtomCoords

@doc """
    MolCharge(mol::GenericMolecule)

Gets the total charge of the molecule.
"""
MolCharge

@doc """
    MolSpin(mol::GenericMolecule)

Gets the total spin of the molecule. Each unpaired electron contributes 1/2.
"""
MolSpin

@doc """
    MolEnergyDataAvailable(mol::GenericMolecule)

Gets the availability of the energy data of the molecule.
"""
MolEnergyDataAvailable

@doc """
    MolEnergyLevels(mol::GenericMolecule [,spin=1|2])

Gets the energy levels of the molecule's MOs.

- `spin`: For closed-shell molecules, the `spin` param should be neglected.
          For open-shell molecules (with non-zero spins), `spin=1` indicates α orbitals and `spin=2` indicates β orbitals, neglecting `spin` would return both two sets of orbitals.
"""
MolEnergyLevels

@doc """
    MolEnergyLevel(mol::GenericMolecule, orbit_ridx)

Gets the energy level of the molecule's selected MO.

- `orbit_ridx`: Index of selected orbit relative to the HOMO (e.g., `0` indicates HOMO, and `-1` indicates HOMO-1).
                For open-shell molecules, according to α/β spins, should be passed in format `(spin, idx)` where for α orbitals `spin=1` and for β orbitals `spin=2`.
"""
MolEnergyLevel

@doc """
    MolOrbitalOccupation(mol::GenericMolecule [,spin=1|2])

Gets the occupation of the molecule's MOs.

- `spin`: For closed-shell molecules, the `spin` param should be neglected.
          For open-shell molecules (with non-zero spins), `spin=1` indicates α orbitals and `spin=2` indicates β orbitals, neglecting `spin` would return both two sets of orbitals.
"""
MolOrbitalOccupation

@doc """
    MolHOMOEnergy(mol::GenericMolecule [,spin=1|2])

Gets the energy of the molecule's HOMO.

- `spin`: For closed-shell molecules, the `spin` param should be neglected.
          For open-shell molecules (with non-zero spins), `spin=1` indicates α orbitals and `spin=2` indicates β orbitals, neglecting `spin` would give both.
"""
MolHOMOEnergy

@doc """
    MolWFATAvailableIndices(mol::GenericMolecule)

Gets the available orbital indices (relative to HOMO) of the molecule's WFAT data.
"""
MolWFATAvailableIndices

@doc """
    MolWFATData(mol::GenericMolecule, orbit_ridx)

Gets the WFAT data in format `(μ, int_data)`.

- `orbit_ridx`: Index of selected orbit relative to the HOMO (e.g., `0` indicates HOMO, and `-1` indicates HOMO-1).
                For open-shell molecules, according to α/β spins, should be passed in format `(spin, idx)` where for α orbitals `spin=1` and for β orbitals `spin=2`.
"""
MolWFATData

@doc """
    MolWFATStructureFactor_G(mol::GenericMolecule, orbit_ridx, nξ, m, θ, χ)

Gets the WFAT structure factor ``G_{n_ξ m}`` according to the given Euler angles `θ` and `χ` (z-y'-z'' convention).
Note: the rotational Euler angles of the molecule would not be applied.

## Parameters
- `orbit_ridx`: Index of selected orbit relative to the HOMO (e.g., `0` indicates HOMO, and `-1` indicates HOMO-1).
                For open-shell molecules, according to α/θ spins, should be passed in format `(spin, idx)` where for α orbitals `spin=1` and for β orbitals `spin=2`.
- `nξ`  : Parabolic quantum number nξ=0,1,2,⋯ (nξ up to 5 is calculated by default).
- `m`   : Parabolic quantum number m=⋯,-1,0,1,⋯ (|m| up to 5 is calculated by default).
- `θ`   : Euler angle θ, passed as a `Real` value or an `AbstractVector` of `Real`.
- `χ`   : Euler angle χ, passed as a `Real` value or an `AbstractVector` of `Real`.
"""
MolWFATStructureFactor_G

@doc """
    MolWFATMaxChannels(mol::GenericMolecule, orbit_ridx)

Gets the maximum value of nξ and |m| calculated in the WFAT integral data.

- `orbit_ridx`: Index of selected orbit relative to the HOMO (e.g., `0` indicates HOMO, and `-1` indicates HOMO-1).
                For open-shell molecules, according to α/β spins, should be passed in format `(spin, idx)` where for α orbitals `spin=1` and for β orbitals `spin=2`.
"""
MolWFATMaxChannels

@doc """
    MolAsympCoeffAvailableIndices(mol::GenericMolecule)

Gets the available orbital indices (relative to HOMO) of the molecule's asymptotic coefficients.
"""
MolAsympCoeffAvailableIndices

@doc """
    MolAsympCoeff(mol::GenericMolecule, orbit_ridx)

Gets the asymptotic coefficients of the molecule.

- `orbit_ridx`: Index of selected orbit relative to the HOMO (e.g., `0` indicates HOMO, and `-1` indicates HOMO-1).
                For open-shell molecules, according to α/β spins, should be passed in format `(spin, idx)` where for α orbitals `spin=1` and for β orbitals `spin=2`.
"""
MolAsympCoeff

@doc """
    MolAsympCoeff_lMax(mol::GenericMolecule, orbit_ridx)

Gets the maximum value of l calculated in the asymptotic coefficients.

- `orbit_ridx`: Index of selected orbit relative to the HOMO (e.g., `0` indicates HOMO, and `-1` indicates HOMO-1).
                For open-shell molecules, according to α/β spins, should be passed in format `(spin, idx)` where for α orbitals `spin=1` and for β orbitals `spin=2`.
"""
MolAsympCoeff_lMax

@doc """
    MolRotation(mol::GenericMolecule) -> (α,β,γ)

Gets the Euler angles (z-y'-z'' convention) specifying the molecule's orientation in format (α,β,γ) (in radian).
"""
MolRotation

@doc """
    SetMolRotation!(mol::GenericMolecule, α,β,γ)
    SetMolRotation!(mol::GenericMolecule, (α,β,γ))

Sets the Euler angles (z-y'-z'' convention) specifying the molecule's orientation in format (α,β,γ) (numerically in radian or a `Unitful.Quantity`).
"""
SetMolRotation!

@doc """
    MolExportAtomInfo(mol::GenericMolecule)

Exports the given molecule's atom information to string as `MolecularCalculator`'s input.
"""
MolExportAtomInfo

@doc """
    MolInitCalculator!(mol::GenericMolecule, MCType::Type=PySCFMolecularCalculator [;kwargs...])

Initializes the `MolecularCalculator` of `mol` with given parameters.
- `MCType`      : Type of `MolecularCalculator` if it is not initialized (default is `PySCFMolecularCalculator`).
- `kwargs...`   : Keyword arguments to pass to the initializer of `MolecularCalculator` of `MCType`, e.g., `basis`, ...
"""
MolInitCalculator!

@doc """
    MolCalcWFATData!(mol::GenericMolecule [,orbit_ridx=0] [;kwargs...])

Calculates the WFAT data of the molecule.
- `orbit_ridx` : Index of selected orbit relative to the HOMO (e.g., `0` indicates HOMO, and `-1` indicates HOMO-1).
                 For open-shell molecules, according to α/β spins, should be passed in format `(spin, idx)` where for α orbitals `spin=1` and for β orbitals `spin=2`.
- `kwargs...` : Keyword arguments to pass to the [`calc_WFAT_data`](@ref) method, e.g. `grid_rNum`, `grid_rMax`, `sf_lMax`, ⋯
"""
MolCalcWFATData!

@doc """
    MolCalcAsympCoeff!(mol::GenericMolecule, orbit_ridx; kwargs...)

Calculates the asymptotic coefficients of the molecule.
- `orbit_ridx` : Index of selected orbit relative to the HOMO (e.g., `0` indicates HOMO, and `-1` indicates HOMO-1).
                 For open-shell molecules, according to α/β spins, should be passed in format `(spin, idx)` where for α orbitals `spin=1` and for β orbitals `spin=2`.
- `kwargs...` : Keyword arguments to pass to the [`calc_asymp_coeff`](@ref), e.g. `grid_rNum`, `l_max`.
"""
MolCalcAsympCoeff!

@doc """
    MolSaveDataAs!(mol::GenericMolecule, data_path [,overwrite=false])

Saves the data of the `GenericMolecule` to the `data_path`. To overwrite the existing file, set `overwrite=true`.
"""
MolSaveDataAs!
