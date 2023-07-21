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

- `laser::Laser`

A `Lasers.Laser` object containing information of the laser field.
Cf. the documentation for [Lasers](@ref lasers_doc).

### Targets

- `target::Target`

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
However, for special demands, the user may refer to the documentation of [`MolCalcMOADKCoeff!`](@ref) and [`Targets.MolecularCalculators.calcMOADKCoeff`](@ref) for more calculation parameters.

#### Calculate WFAT Data

To calculate the data necessary for the [WFAT](@ref WFAT) of the `Molecule` object, invoke the [`MolCalcWFATData!`](@ref) method:
```julia
MolCalcWFATData!(mol, orbitIdx_relHOMO = 0)
```

To obtain the data of other orbitals besides HOMO, the user may alter the `orbitIdx_relHOMO` parameter.
For custom calculation parameters, refer to the documentation of [`MolCalcWFATData!`](@ref) and [`Targets.MolecularCalculators.calcStructFactorData`](@ref).

#### Setting the Molecule's Orientation

The orientation of the molecule can be specified by invoking the [`SetMolRotation`](@ref) method:
```julia
SetMolRotation(mol, 0.0,π/2,π/3)
```
and can be obtained through the [`MolRotation`](@ref) method:
```julia
MolRotation(mol)
```



## Initial Condition Methods

- `init_cond_method = <:ADK|:SFA|:SFAAE|:WFAT|:MOADK>`

Method of electrons' initial conditions.
Currently supports `:ADK`, `:SFA`, `:SFAAE` for atoms, and `:WFAT`, `:MOADK` for molecules.

For more information about the theories, cf. [Theory - Initial Conditions](@ref theory_init_cond).

!!! note "Note"
    For initial condition methods [ADK](@ref ADK), [SFA](@ref SFA) or [SFA-AE](@ref SFAAE), you must specify an atom target of types [`HydrogenLikeAtom`](@ref) or [`SAEAtom`](@ref);
    For initial condition methods [WFAT](@ref WFAT) or [MO-ADK](@ref MOADK), you must specify a molecule target of type [`Molecule`](@ref).

    |  | [ADK](@ref ADK) | [SFA](@ref SFA) | [SFA-AE](@ref SFAAE) | [WFAT](@ref WFAT) | [MO-ADK](@ref MOADK) |
    | :------------------------- |:-:|:-:|:-:|:-:|:-:|
    | [`HydrogenLikeAtom`](@ref) | ✔ | ✔ | ✔ |   |   |
    | [`SAEAtom`](@ref)          | ✔ | ✔ | ✔ |   |   |
    | [`Molecule`](@ref)         |   |   |   | ✔ | ✔ |

### Atomic Rate Prefix

- `rate_prefix = <:ExpRate|:ExpPre|:ExpJac|:Full>`

Prefix of the exponential term in the ionization rate (default `:ExpRate`).

For atomic [ADK](@ref ADK), [SFA](@ref SFA) and [SFA-AE](@ref SFAAE), we obtained the ionization probability in the following form:
```math
\mathrm{d}W/\mathrm{d}\bm{k}_\perp \mathrm{d}t_{\mathrm{r}} = \bm{J}(k_d,t_{\mathrm{r}}) \lvert P_{\bm{p}}(t_{\mathrm{s}}) \rvert^2 \exp(2\ \mathrm{Im}\ \Phi_{\mathrm{tun}}),
```
where the ``\lvert P_{\bm{p}}(t_{\mathrm{s}}) \rvert^2`` denotes the prefactor, which has different expressions for different theories:

| [SFA](@ref SFA) | [SFA-AE](@ref SFAAE) | [ADK](@ref ADK) |
| :-: | :-: | :-: |
| ``\{ [\bm{p}+\bm{A}(t_{\mathrm{s}})] \cdot \bm{F}(t_{\mathrm{s}}) \}^{-\alpha}`` | ``\left[ (k_\perp^2+2I_{\mathrm{p}})(F^2-\bm{k}_\perp \cdot \bm{F}') \right]^{-\alpha/2}`` | ``\left[ (k_\perp^2+2I_{\mathrm{p}})F^2\right]^{-\alpha/2}`` |

The ``\bm{J}(k_d,t_{\mathrm{r}})`` denotes the Jacobian which arises from the coordinate transformation.

For initial condition methods [ADK](@ref ADK), [SFA](@ref SFA) and [SFA-AE](@ref SFAAE),
specifying `ExpRate` would not add any prefix besides the exponential term;
`ExpPre` would include the prefactor ``\lvert P_{\bm{p}}(t_{\mathrm{s}}) \rvert^2``;
`ExpJac` would include the Jacobian ``\bm{J}(k_d,t_{\mathrm{r}})``;
`Full` indicates inclusion of both the prefactor and Jacobian.

### Tunneling Exit Methods For Atomic ADK

- `adk_tun_exit = <:IpF|:FDM|:Para>`

Tunneling exit method for atomic ADK methods (when `init_cond_method==:ADK`) (default `:IpF`).
Cf. [Tunneling Exit Methods For Atomic ADK](@ref tun_exit_atomic_adk).



## Sampling Methods and Parameters

In *SemiclassicalSFI.jl* the initial electrons are sampled in the ``(t_{\mathrm{r}},k_d,k_z)`` coordinate, where ``k_d`` denotes the initial momentum's component in the ``xy`` plane (which is perpendicular to the electric field).
There are two ways of sampling in these coordinates, namely *step sampling* and *Monte-Carlo sampling*.

- `sample_monte_carlo = false` : Determines whether Monte-Carlo sampling is used when generating electron samples (default `false`).

- `sample_t_interval = (start,stop)` : Time interval in which the initial electrons are sampled.

- `sample_t_num` : Number of time samples.

### Step Sampling

In the step sampling scheme, the `sample_t_num` time samples are uniformly distributed in the interval `sample_t_interval`.
In each time sample, a batch of electrons of different initial conditions are launched and collected, whose initial momenta ``\bm{k}_\perp`` are distributed on a Cartesian grid ``(k_d,k_z)``.
The Cartesian grid is defined by the following parameters:

- `ss_kd_max`, `ss_kd_num`, `ss_kz_max`, `ss_kz_num`

In the ``k_d`` dimension, `ss_kd_num` samples distribute uniformly in the interval (-`ss_kd_max`,`ss_kd_max`);
and in the ``k_z`` dimension, there are `ss_kz_num` equidistant samples in the interval (-`ss_kz_max`,`ss_kz_max`).

The step sampling method is supported for all initial condition methods.

### Monte-Carlo Sampling

In the Monte-Carlo sampling scheme, the `sample_t_num` time samples are randomly chosen in the `sample_t_interval`;
A batch containing `mc_t_batch_size` electrons would be sampled in a single time sample, the electrons' initial momenta ``\bm{k}_\perp`` are also randomly sampled inside a circle ``k_d^2+k_z^2 \leq k_{\perp\mathrm{max}}^2``, where the ``k_{\perp\mathrm{max}}^2`` is defined in the parameter as `mc_kt_max`.

- `mc_t_batch_size` : Number of electron samples in a single time sample.
- `mc_kt_max` : Maximum value of momentum's transversal component (perpendicular to field direction).

## Trajectory Simulation

### Phase Methods

### Non-dipole Effect

### Accuracy Control

### GPU Acceleration (Experimental)


## Final Electron Collecting & Saving

### 2D/3D Momentum Spectrum Collecting

### Rydberg Final State Collecting

### Output File Name