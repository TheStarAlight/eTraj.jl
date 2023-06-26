# Targets

*This section provides information of available targets (atoms/molecules) in the library.*

A target interacts with the laser field and release the electron through multi-photon or tunneling processes.
Here we list available targets implemented in the [`Targets`](@ref) module of the library.

```@docs
Targets
```

```@contents
Pages = ["manual1_targets.md"]
Depth = 3
```

```@meta
CurrentModule = SemiclassicalSFI.Targets
```

## Hydrogen-Like Atom

A hydrogen-like atom has a potential of the form
```math
V(r) = -\frac{Z}{\sqrt{r^2+a}},
```
where ``Z`` is the nuclear charge number;
``a`` denotes the soft-core parameter, which is applied to avoid singularity of the potential and can be adjusted to fit the actual ionization potential of the atom (obtained by TDSE).

The hydrogen-like atom is implemented in the library as [`HydrogenLikeAtom`](@ref).
```@docs
Targets.HydrogenLikeAtom
```


## Single-Active-Electron (SAE) Atom

The single-active-electron (SAE) atom is an implementation of the empirical atomic SAE model potential proposed by Tong *et al.* [^Tong_2005]
The model potential of the SAE atom has the form
```math
V(r) = -\frac{Z + a_1 \mathrm{e}^{-b_1 r} + a_2 r \mathrm{e}^{-b_2 r} + a_3 \mathrm{e}^{-b_3 r}}{r},
```
where the ``a_i`` and ``b_i`` are tunable model potential parameters [^note].

The SAE atom is implemented in the library as [`SAEAtom`](@ref).
```@docs
Targets.SAEAtom
```

[^Tong_2005]: X. M. Tong *et al.*, Empirical Formula for Static Field Ionization Rates of Atoms and Molecules by Lasers in the Barrier-Suppression Regime. *J. Phys. B: At. Mol. Opt. Phys.* **38**, 2593–2600. DOI: [10.1088/0953-4075/38/15/001](https://dx.doi.org/10.1088/0953-4075/38/15/001)
[^note]: The symbols of the parameters are different from that in the original article. ``a_1, b_1, a_2, b_2, a_3, b_3`` correspond to ``a_1, a_2, a_3, a_4, a_5, a_6`` in the original article respectively.


## Preset Atom Library

The library provides some preset commonly-used atoms or atomic ions for convenience.

- [`HydrogenLikeAtom`](@ref) : H, He⁺, Li²⁺

- [`SAEAtom`](@ref) :          He, Ne, Ne⁺, Ne²⁺, Ar, Ar⁺, Ar²⁺, V, Ni, Kr, Kr⁺, Rb, Nb, Pd, Xe, Xe⁺, Ta

These atoms/ions can be obtained by invoking `Targets.**Atom()` (for neutral atoms) or `Targets.**#pAtom()` (for positive atomic ions), where `**` denotes the symbol of the nucleus and `#` denotes the positive charge the ion carries.

Example:

```@setup manual_targets
using SemiclassicalSFI
```
```@repl manual_targets
t1 = Targets.HAtom()
t2 = Targets.Xe1pAtom()
```


## Molecule

The `Molecule` object represents a generic molecule, which is implemented in the library as [`Molecule`](@ref).
The structure of `Molecule` is much more complex than that of atoms because the [Molecular ADK (MO-ADK)](@ref) and [Weak-Field Asymptotic Theory (WFAT)](@ref) features for molecular strong-field ionization require a number of coefficients, which are saved to files for convenience.

### Initialization, saving and loading

The `Molecule` object can be initialized either by providing necessary information of the molecule (mainly atoms, coordinates of the atoms and the charge of the molecule) or from external data (stored in the HDF5 format), cf. the documentation of [`Molecule`](@ref):

```@docs
Targets.Molecule
```

The molecule object, after modification, can be manually saved to a HDF5 file via [`MolSaveDataAs`](@ref).

```@docs
Targets.MolSaveDataAs
```

### Molecular-SFI Data Preparation

To use the molecular strong-field ionization theories such as the MO-ADK and WFAT to provide the intitial conditions of the electrons, the structure coefficients of the `Molecule` have to be calculated beforehand and stored in the object.
Cf. the documentation of [`MolCalcMOADKCoeff!`](@ref) and [`MolCalcWFATData!`](@ref).

Evaluation of the structure coefficients depends on the external quantum chemistry packages.
The [`Targets.MolecularCalculators`](@ref) module undertakes the task of communication with the external quantum packages.
Currently only the [`PySCFMolecularCalculator`](@ref) is implemented.

!!! note "Customized calculation parameters"
    When invoking `MolCalcMOADKCoeff!` and `MolCalcWFATData!` to perform calculation of structure coefficients, customized calculation parameters can be passed to the `kwargs` of these methods.
    These parameters would be passed to the constructor method of the `MolecularCalculator` (e.g., the `basis` parameter of the [`PySCFMolecularCalculator`](@ref)), as well as the [`MolecularCalculators.calcStructFactorData`](@ref), [`MolecularCalculators.calcMOADKCoeff`](@ref) methods.
    Refer to their documentation below for more information.

```@docs
Targets.MolCalcMOADKCoeff!
Targets.MolCalcWFATData!
```

```@docs
Targets.MolecularCalculators
Targets.MolecularCalculators.PySCFMolecularCalculator
```

```@docs
Targets.MolecularCalculators.calcStructFactorData
Targets.MolecularCalculators.calcMOADKCoeff
```

### Molecule's Orientation

The molecule's orientation is described by a set of Euler angles (``z-y'-z''`` convention), which defines a rotational transformation from the molecular frame (MF) to the lab frame (LF).
This property of `Molecule` is NOT included in the saved file and thus needs to be specified each time upon initialization of the `Molecule` object from external files.

The orientation of the molecule can be obtained and set via the [`MolRotation`](@ref) and [`SetMolRotation`](@ref) methods.

```@docs
Targets.MolRotation
Targets.SetMolRotation
```
