# [Manual: The `Targets` Module](@id targets_doc)

The [`Targets`](@ref) module implements the abstraction of targets and provides parameters of some commonly-used targets.

```@docs
Targets
Target
```

```@contents
Pages = ["manual2_targets.md"]
Depth = 3
```

---------------------------

## Atoms under SAE approximation

The [`HydrogenLikeAtom`](@ref) and [`SAEAtom`](@ref) are both subtypes of the [`SAEAtomBase`](@ref) base type, which represents an atom under the SAE approximation.
A key ingredient of atomic objects lies in the potential function of the residual ion after the electron gets ionized, which is the only difference between the two types.

The [`HydrogenLikeAtom`](@ref)'s potential function is of the form:
```math
\begin{equation}
    V(r) = -\frac{Z}{\sqrt{r^2+a}},
\end{equation}
```
with $a$ the soft-core parameter which avoids singularity of the potential during numerical simulation.

The [`SAEAtom`](@ref)'s potential function is adopted from Tong's model [^Tong_2005]:
```math
\begin{equation}
    V(r) = -\frac{Z+a_1\ee^{-b_1 r}+a_2 r\ee^{-b_2 r}+a_3\ee^{-b_3 r}}{\sqrt{r^2+a}},
\end{equation}
```
where $a_i$ and $b_i$ are tunable parameters to fit the effective potential felt by the electron.

[^Tong_2005]: X. M. Tong and C. D. Lin, Empirical formula for static field ionization rates of atoms and molecules by lasers in the barrier-suppression regime, *J. Phys. B: At. Mol. Opt. Phys.* **38**, 2593 (2005). DOI: [10.1088/0953-4075/38/15/001](https://doi.org/10.1088/0953-4075/38/15/001)

```@docs
SAEAtomBase
HydrogenLikeAtom
HydrogenLikeAtom()
SAEAtom
SAEAtom()
```

For the convenience of the user, there are some presets of commonly-used atoms, which can be accessed through the [`get_atom`](@ref) method.
Available keys are accessed by invoking [`get_available_atoms`](@ref).

- [`HydrogenLikeAtom`](@ref) : H, He⁺, Li²⁺
- [`SAEAtom`](@ref) :          He, Ne, Ne⁺, Ne²⁺, Ar, Ar⁺, Ar²⁺, V, Ni, Kr, Kr⁺, Rb, Nb, Pd, Xe, Xe⁺, Ta

```@docs
get_atom
get_available_atoms
```

---------------------------

## Molecules

For molecule targets, the [`GenericMolecule`](@ref) type is implemented based on the [`MoleculeBase`](@ref) supertype, whose structure is more complicated.
A [`GenericMolecule`](@ref) stores information about the atoms that make up the molecule, together with their coordinates,
as well as the asymptotic coefficients [``C_{lm}`` in Eq. (38) in [Molecular SFA](@ref MOSFA)] and WFAT's integral coefficients [Eq. (56) in [WFAT](@ref WFAT)],
which are obtained using other quantum chemistry packages.

There are two ways to initialize a [`GenericMolecule`](@ref): build from zero ([`GenericMolecule()`](@ref)) or from an existing file ([`LoadMolecule`](@ref)).

```@docs
MoleculeBase
GenericMolecule
GenericMolecule()
LoadMolecule
```

We also provide some pre-defined molecules, which can be accessed through the [`get_mol`](@ref) method.

```@docs
get_mol
get_available_mols
```

### Molecule's Orientation

The molecule's orientation is described by a set of Euler angles (``z-y'-z''`` convention), which defines a rotational transformation from the molecular frame (MF) to the lab frame (LF).
This property is NOT included in the saved file and thus needs to be specified each time upon initialization of the [`GenericMolecule`](@ref) object from external files.

!!! note "Note"
    Here the three Euler angles `(α,β,γ)` that describe the [`GenericMolecule`](@ref)'s orientation are completely different from that of the Euler angles `(θ,χ)` in the [WFAT](@ref WFAT) and [MO-SFA](@ref MOSFA) theory.
    These theories' "lab frame" is chosen for convenience of theoretical formulation, where the electric field is assumed to be static, pointing towards the ``+z`` direction,
    and has no relation with the lab frame mentioned above.

The orientation of the molecule can be obtained and set via the [`MolRotation`](@ref) and [`SetMolRotation!`](@ref) methods.

### Quantum Chemistry Calculation

As for the quantum chemistry calculation, we implemented the scheme in `PySCFMolecularCalculator` using the [`PySCF`](https://github.com/pyscf/pyscf), which works on Linux and macOS platforms (as a remedy, Windows users can use the [Windows Subsystem of Linux (WSL)](https://github.com/microsoft/WSL)).
The calculation scheme of WFAT structure factor is adopted from [`PyStructureFactor`](https://github.com/TheStarAlight/PyStructureFactor).
Future extension is possible by implementing the supertype `MolecularCalculatorBase`.
Since there are some presets of molecules available via [`get_mol`](@ref), we are not going to detail on the manual of running the calculation in the text.

The example of initializing and calculating the essential data of [`GenericMolecule`](@ref) in REPL is presented as follows:
```julia-repl
julia> using eTraj.Targets, eTraj.Units

julia> mol = GenericMolecule(atoms=["O","C","O"], atom_coords=[0 0 -1.1600; 0 0 0; 0 0 1.1600]*Å, charge=0, name="Carbon Dioxide (CO₂)")
[GenericMolecule] Carbon Dioxide (CO₂)

julia> MolInitCalculator!(mol, basis="cc-pVTZ")
[ Info: [PySCFMolecularCalculator] Running molecular calculation...

julia> MolCalcAsympCoeff!(mol, 0); MolCalcAsympCoeff!(mol, -1)
[ Info: [PySCFMolecularCalculator] Running calculation of asymptotic coefficients... (ionizing orbital HOMO)
[ Info: [PySCFMolecularCalculator] Running calculation of asymptotic coefficients... (ionizing orbital HOMO-1)

julia> MolCalcWFATData!(mol, 0); MolCalcWFATData!(mol, -1)
[ Info: [PySCFMolecularCalculator] Running calculation of WFAT structure factor data... (ionizing orbital HOMO)
[ Info: [PySCFMolecularCalculator] Running calculation of WFAT structure factor data... (ionizing orbital HOMO-1)

julia> MolSaveDataAs!(mol, "Molecule_CO2.jld2")
[ Info: [GenericMolecule] Data saved for molecule Carbon Dioxide (CO₂) at `Molecule_CO2.jld2`.

julia> mol_ = LoadMolecule("Molecule_CO2.jld2")  # load from saved file
[GenericMolecule] Carbon Dioxide (CO₂)
Asymp coeff of HOMO-1 & HOMO available
WFAT data of HOMO-1 & HOMO available
#          E (Ha)  occp
⋮    ⋮       ⋮      ⋮⋮
13 LUMO+1   0.207  ----
12 LUMO     0.175  ----
11 HOMO    -0.542  -↿⇂-
10 HOMO-1  -0.542  -↿⇂-
9  HOMO-2  -0.714  -↿⇂-
8  HOMO-3  -0.714  -↿⇂-
⋮    ⋮       ⋮      ⋮⋮
```

---------------------------

## List of Available Properties & Methods

### Generic `Target`

```@docs
IonPotential
AsympNuclCharge
TargetName
TargetPotential
TargetForce
TrajectoryFunction
```

---------------------------

### `SAEAtomBase`

```@docs
AngularQuantumNumber
MagneticQuantumNumber
AsympCoeff
SoftCore
QuantizationAxisOrientaion
```

---------------------------

### `GenericMolecule`

```@docs
MolAtoms
MolAtomCoords
MolCharge
MolSpin
MolEnergyDataAvailable
MolEnergyLevels
MolEnergyLevel
MolOrbitalOccupation
MolHOMOEnergy
MolWFATAvailableIndices
MolWFATData
MolWFATStructureFactor_G
MolWFATMaxChannels
MolAsympCoeffAvailableIndices
MolAsympCoeff
MolAsympCoeff_lMax
MolRotation
SetMolRotation!
MolExportAtomInfo
MolInitCalculator!
MolCalcWFATData!
MolCalcAsympCoeff!
MolSaveDataAs!
```

---------------------------

### `MolecularCalculators`

```@docs
MolecularCalculatorBase
PySCFMolecularCalculator
PySCFMolecularCalculator()
Targets.calc_WFAT_data
Targets.calc_asymp_coeff
```
