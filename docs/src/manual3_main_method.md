# Main Method `performSFI`

*This section introduces the main method* [`performSFI`](@ref)*, which performs the trajectory simulation and saves the final electron momentum spectrum.*

```@contents
Pages = ["manual3_main_method.md"]
Depth = 4
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

### Workflow for the `Molecule` Target



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
