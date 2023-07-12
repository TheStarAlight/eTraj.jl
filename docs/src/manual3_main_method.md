# Main Method `performSFI`

*This section introduces the main method* [`performSFI`](@ref)*, which performs the trajectory simulation and saves the final electron momentum spectrum.*

```@contents
Pages = ["manual3_main_method.md"]
Depth = 3
```

```@meta
CurrentModule = SemiclassicalSFI
```


## Brief Documentation

```@docs
performSFI
```

## Lasers & Targets

### Lasers

- `laser:Laser`

A `Lasers.Laser` object containing information of the laser field.
Cf. the documentation for [Lasers](@ref lasers_doc).

### Targets

- `target:Target`

A `Targets.Target` object containing information of the atom/molecule target.
Cf. the documentation for [Targets](@ref targets_doc).

### Workflow for Preparation of the `Molecule` Target

To use the [`Molecule`](@ref) as the target and the relating initial condition methods [MO-ADK](@ref MOADK) or [WFAT](@ref WFAT), some coefficients need to be calculated in advance.
Here we present a workflow for preparation of the `Molecule` target before invoking the `performSFI` method, taking the carbon monoxide molecule as an example.

!!! note "Note"
    Currently only the PySCF is implemented as the library's molecular computation interface, which only supports the Linux platform.
    What's more, in current version, only *close-shell* molecules are supported!


#### Initialization

First of all, initialize a `Molecule` object, and provide the necessary information of the molecule. Cf. the documentation of [`Molecule`](@ref).

```julia
using SemiclassicalSFI.Targets
mol = Molecule(atoms=["C","O"], atom_coords=[0 0 -0.180; 0 0 0.950],
               charge=0, name="Carbon Monoxide",
               data_path="./Molecule_CarbonMonoxide.h5")
```

!!! note "Data saving of the Molecule object"
    If the user specifies `data_path` in the constructor method of `Molecule`, the data would be automatically saved each time the user invokes the [`MolCalcMOADKCoeff!`](@ref) and [`MolCalcWFATData!`](@ref).
    However, if doesn't specify (in case the user does not wish to save the data), the data would not be saved, and the user has to manually invoke [`MolSaveDataAs`](@ref) to save the data afterwards.


#### Calculate MO-ADK Coefficients

To calculate the [MO-ADK](@ref MOADK) coefficients of the `Molecule` object, invoke the [`MolCalcMOADKCoeff!`](@ref) method:
```julia
MolCalcMOADKCoeff!(mol)
```

For typical small molecules, using default parameters usually gives satisfactory results.
However, for special demands, the user may refer to the documentation of [`MolCalcMOADKCoeff!`](@ref) and [`Targets.MolecularCalculators.calcMOADKCoeff`](@ref) for more configuration parameters.

#### Calculate WFAT Data

To calculate the data necessary for the [WFAT](@ref WFAT) of the `Molecule` object, invoke the [`MolCalcWFATData!`](@ref) method:
```julia
MolCalcWFATData!(mol, orbitIdx_relHOMO = 0)
```

To obtain the data of other orbitals besides HOMO, the user may alter the `orbitIdx_relHOMO` parameter.
If the user requires custom calculation parameters, refer to the documentation of [`MolCalcWFATData!`](@ref) and [`Targets.MolecularCalculators.calcStructFactorData`](@ref).


## Initial Condition Methods

### Atomic SFA, SFA-AE and ADK

#### Rate Prefactor

#### Tunneling Exit Methods For ADK


### Molecular WFAT and MO-ADK


## Sampling Methods and Parameters

### Step-Sampling

### Monte-Carlo-Sampling


## Trajectory Simulation

### Phase Methods

### Non-dipole Effect

### Accuracy Control

### GPU Acceleration


## Final Electron Collecting & Saving

### 2D/3D Momentum Spectrum Collecting

### Rydberg Final State Collecting

### Output File Name
