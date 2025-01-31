# Step-by-step Tutorial

*-- For beginners in julia*

This section provides a step-by-step tutorial on how to install and use *eTraj*.

```@contents
Pages = ["stepbystep_tutorial.md"]
Depth = 4
```

-----------------------

## Installing Julia

The manual for installation of Julia is available at the [official website](https://julialang.org/downloads/).

### Windows

For Windows users, install using the following command:
```bash
> winget install julia -s msstore
```
This would also install [`juliaup`](https://github.com/JuliaLang/juliaup) which helps in managing multiple versions of Julia.

The [`jill`](https://github.com/abelsiqueira/jill) is an alternative to `juliaup` which can be used for installing and managing multiple versions of Julia.
It would automatically search for the fastest mirror and download the installer and is useful for users with slow internet connections.

For offline installation, an installer is available at [here](https://julialang.org/downloads/#manual_download).
Please follow the [instructions](https://julialang.org/downloads/platform/#windows) and make sure that julia is added to the `PATH` environment variable.

### Linux & macOS

For Linux & macOS users, use the following command:
```bash
$ curl -fsSL https://install.julialang.org | sh
```
which would install [`juliaup`](https://github.com/JuliaLang/juliaup) for multiple version management.

The [`jill`](https://github.com/abelsiqueira/jill) is an alternative to `juliaup` which can be used for installing and managing multiple versions of Julia.
It would automatically search for the fastest mirror and download the installer and is useful for users with slow internet connections.

For offline installation, the binaries are available at [here](https://julialang.org/downloads/#manual_download).
Please follow the instructions for [Linux](https://julialang.org/downloads/platform/#linux_and_freebsd) and [macOS](https://julialang.org/downloads/platform/#macos) and make sure that julia is added to the `PATH` environment variable.

!!! warning "Note"
    It is not recommended to install julia using Linux package managers or macOS homebrew.
    See this [notice](https://julialang.org/downloads/#please_do_not_use_the_version_of_julia_shipped_by_linux_or_bsd_package_managers).

-----------------------

## Package Management

Configuration files and packages are stored in `~/.julia/` directory.

The package manager is accessible through the REPL[^1] by typing `]`.
This will open a new prompt (starting by `(@v1.x) pkg>`) where you can install (`add`), update (`up`), or remove (`rm`) packages.
To install/update/remove a package, simply type `add`/`up`/`rm` in the package mode, following the package name, and press `Enter`.
For example:
```
julia> ]
(@v1.9) pkg> add Example
   Resolving package versions...
   Installed Example ─ v0.5.3
    Updating `~/.julia/environments/v1.9/Project.toml`
  [7876af07] + Example v0.5.3
    Updating `~/.julia/environments/v1.9/Manifest.toml`
  [7876af07] + Example v0.5.3
```
To view the list of installed packages, type `status` or `st` in package mode and press `Enter`.

See the [documentation](https://docs.julialang.org/en/v1/stdlib/Pkg/) for more details on using the package manager.

[^1]: REPL (Read–eval–print loop) is an interactive programming environment that reads user input, evaluates it, and prints the result. The Julia REPL is opened by running `julia` in your terminal.

To install `eTraj`, use the following command in package mode:
```
(@v1.9) pkg> add https://github.com/TheStarAlight/eTraj.jl.git
```

For offline installation, use git to clone the repository:
```
$ git clone https://github.com/TheStarAlight/eTraj.jl.git
```
and then install it in Julia by specifying the directory:
```
(@v1.9) pkg> add /path/to/eTraj
```

-----------------------

## [Configuring Python and PySCF (*Optional*)](@id config_py)

Currently, the calculation of molecules' asymptotic coefficients (for MO-ADK/MO-SFA) and WFAT coefficients rely on the [`PySCF`](https://github.com/pyscf/pyscf) python package. `eTraj` calls the PySCF using the [`PyCall.jl`](https://github.com/JuliaPy/PyCall.jl) package.

There are two ways to set up the Python environment used by `PyCall`:

1. using your local Python environment by specifying the path of your Python executable in `ENV["PYTHON"]` and build the PyCall package.
2. using a private Python environment managed by the [`Conda.jl`](https://github.com/JuliaPy/Conda.jl), which is implicitly installed by the `PyCall` package by default;

### Using the local Python environment

To correctly set up the configuration of `PyCall`, first, set the `PYTHON` environment variable to the path your Python executable, and build the `PyCall` package:

```julia
ENV["PYTHON"] = "path/to/python_exec"
using Pkg
Pkg.build("PyCall")
# In pkg mode of REPL:
# (@v1.9) pkg> build PyCall
```

And don't forget to install `PySCF` via pip in your system shell:

```bash
$ pip install pyscf==2.3.0
```

### Using Conda.jl

Before installing `eTraj`, install `Conda` first:

```julia
using Pkg
Pkg.add("Conda")
```

Then call `pip` within `Conda` to install `PySCF`:

```julia
using Conda
Conda.pip_interop(true)
Conda.pip("install", "pyscf==2.3.0")
```

!!! warning "Note"
    Since the `PySCF` does not support the Windows platform, the molecular calculation must be performed on a Linux or macOS platform.
    However, for Windows users, they may install the [WSL](https://learn.microsoft.com/en-us/windows/wsl/install) (Windows Subsystem for Linux), which supports the `PySCF`.

-----------------------

## [Running Scripts](@id running-scripts)

To run a julia script, use the following command in your terminal:
```bash
$ julia /path/to/script.jl
```
To enable multithreading, use the `-t` option followed by the number of threads you want to use.
For example, the following command runs a script with 4 threads:
```bash
$ julia -t 4 /path/to/script.jl
```

If you want to run a script interactively, you can start Julia in REPL mode and then load your script using the `include` function.
For example:
```julia
julia> include("/path/to/script.jl")
```

-----------------------

## Using `eTraj`

It is suggested to write the simulation task in a julia script, and then run it using julia.

!!! note "Tips"
    It is recommended to use [vscode](https://code.visualstudio.com/) with [julia extension](https://github.com/julia-vscode/julia-vscode) installed for GUI or the vim editor with [julia-vim](https://github.com/JuliaEditorSupport/julia-vim) installed for TUI in order for a better coding experience.
    The julia extensions provide the feature of typing special Unicode characters in a LaTeX-like syntax,
    e.g., if you want to type the Greek letter α, you can type `\alpha`, and press `Tab(⇥)`.
    This feature is also available in Julia REPL.

-----------------------

### Defining Lasers

The following example initializes a 4-cycle linearly polarized 800-nm laser with its peak intensity at 0.1 PW/cm² (1.0×10¹⁴ W/cm²) (See also [`Cos2Laser`](@ref) & [`Cos4Laser`](@ref)):
```julia
using eTraj.Lasers, eTraj.Units
l = Cos2Laser(                # Cos²-shape in field strength's temporal profile
    peak_int = 1e14W/cm^2,    # Peak intensity (`0.1PW/cm^2` also works)
    wave_len = 800.0nm,       # wavelength at 800 nm (`0.8μm` also works)
    cyc_num  = 4,             # 10 cycles in duration (possible substitution: `duration = 10.6fs`)
    cep      = 0.0,           # carrier-envelope phase
    ellip    = 0.0            # linear polarization
)
```
Another example for Gaussian-shaped laser (See also [`GaussianLaser`](@ref)):
```julia
using eTraj.Lasers, eTraj.Units
l = GaussianLaser(
    peak_int        = 1.0PW/cm^2,
    ang_freq        = 40.0eV,
    FWHM_duration   = 10.0fs,
    ellip           = 1.0
)
```
It is also possible to combine two monochromatic lasers to form a [`BichromaticLaser`](@ref):
```julia
using eTraj.Lasers, eTraj.Units
l1 = Cos4Laser(peak_int=1.0PW/cm^2, wave_len=800nm, cyc_num=10, ellip=1)
l2 = Cos4Laser(peak_int=1.0PW/cm^2, wave_len=400nm, cyc_num=20, ellip=-1)
l = BichromaticLaser(l1=l1, l2=l2, delay=0.5fs)
```

-----------------------

### Defining Targets

`eTraj` provides flexible interfaces for defining targets.

#### Atoms

For atom targets, you can use either predefined atoms or create custom ones.

A [`HydrogenLikeAtom`](@ref)'s potential is a soft-core Coulomb potential.
An [`SAEAtom`](@ref)'s potential is defined according to Tong's model [[*JPB* **38**, 2593 (2005)](https://doi.org/10.1088/0953-4075/38/15/001)].
An example of initialization:
```julia
using eTraj.Targets, eTraj.Units
t1 = HydrogenLikeAtom(
    Ip = 0.5,           # ionization potential in atomic units (or `13.6eV`)
    Z  = 1,             # nuclear charge
    soft_core = 1e-10,  # soft-core parameter to avoid singularity
    name = "Custom H"   # optional name
)
t2 = SAEAtom(
    Ip = 12.13eV,
    Z  = 1,
    l  = 1,         # angular momentum quantum number
    a1 = 51.35554,  # model parameters for effective potential
    b1 = 2.111554,
    a2 = -99.92747,
    b2 = 3.737221,
    a3 = 1.644457,
    b3 = 0.4306465,
    asymp_coeff = 1.3,  # asymptotic coefficient of wavefunction
    name = "Custom Xe"  # optional name
)
```

To define a target using predefined atoms, use the [`get_atom`](@ref) method:

```jldoctest
julia> using eTraj.Targets, eTraj.Units

julia> t1 = get_atom("H")      # Hydrogen atom
[HydrogenLikeAtom] Atom H, Ip=0.5000 (13.61 eV), Z=1

julia> t2 = get_atom("He1p")   # He⁺ ion
[HydrogenLikeAtom] Atom He⁺, Ip=1.0000 (27.21 eV), Z=2

julia> t3 = get_atom("Xe"; m=1, quan_ax_θ=π/2, quan_ax_ϕ=0.0) # Xenon atom with magnetic quantum number m=1 (quantization axis is x-axis)
[SAEAtom] Atom Xe (p orbital, m=1), Ip=0.4458 (12.13 eV), Z=1, θϕ=(90.0°,0.0°)
```

You can view available predefined atoms using the [`get_available_atoms`](@ref) method:
```jldoctest
julia> using eTraj.Targets

julia> get_available_atoms() |> print
["H", "He1p", "Li2p", "He", "Ne", "Ne1p", "Ne2p", "Ar", "Ar1p", "Ar2p", "V", "Ni", "Kr", "Kr1p", "Rb", "Nb", "Pd", "Xe", "Xe1p", "Ta"]
```

#### Molecules

Some commonly used molecules are predefined and can be accessed using the [`get_mol`](@ref) method:
```julia-repl
julia> using eTraj.Targets, eTraj.Units

# H₂ molecule with z-axis as the molecular axis
julia> m1 = get_mol("Hydrogen")
[GenericMolecule] Hydrogen (H₂)
Asymp coeff of HOMO available
WFAT data of HOMO available
#          E (Ha)  occp
⋮    ⋮       ⋮      ⋮⋮
3  LUMO+1   0.232  ----
2  LUMO     0.148  ----
1  HOMO    -0.594  -↿⇂-

# O₂ molecule with a non-zero total spin
julia> m2 = get_mol("Oxygen")
[GenericMolecule] Oxygen (O₂)
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

# CO molecule rotated 90° around y axis
julia> m3 = get_mol("Carbon Monoxide", rot_β=90°)
[GenericMolecule] Carbon Monoxide (CO), αβγ=(0.0°,90.0°,0.0°)
Asymp coeff of HOMO-2 & HOMO-1 & HOMO available
WFAT data of HOMO-2 & HOMO-1 & HOMO available
#          E (Ha)  occp
⋮    ⋮       ⋮      ⋮⋮
9  LUMO+1   0.139  ----
8  LUMO     0.139  ----
7  HOMO    -0.554  -↿⇂-
6  HOMO-1  -0.639  -↿⇂-
5  HOMO-2  -0.639  -↿⇂-
4  HOMO-3  -0.804  -↿⇂-
⋮    ⋮       ⋮      ⋮⋮
```

You can view available predefined atoms and molecules using the [`get_available_mols`](@ref) method:
```jldoctest
julia> using eTraj.Targets

julia> get_available_mols() |> print
["Hydrogen", "Nitrogen", "Oxygen", "Carbon Monoxide", "Nitric Oxide", "Hydrochloric Acid", "Carbon Dioxide", "Sulfur Dioxide", "Water", "Ammonia", "Acetylene", "Methane", "Benzene"]
```

To create a custom molecule, use the [`GenericMolecule`](@ref) constructor.
The following example creates a carbon dioxide molecule, and runs quantum chemistry calculations to obtain the asymptotic coefficients (for MO-SFA-based initial condition methods) and WFAT coefficients:
```julia
using eTraj.Targets, eTraj.Units

# Creating a custom molecule
m = GenericMolecule(
    atoms = ["O", "C", "O"],    # atom types
    atom_coords = [             # atomic coordinates
        0 0 -1.1600;
        0 0  0.0000;
        0 0  1.1600
    ] * Å,
    charge = 0,                 # total charge
    name = "Carbon Dioxide"     # optional name
)

# Calculating asymptotic coefficients and WFAT data
MolInitCalculator!(m, basis="cc-pVDZ")  # initializes the calculator backend with specified basis set
MolCalcAsympCoeff!(m, 0)   # calculates asymptotic coefficients for HOMO (index 0)
MolCalcWFATData!(m, 0)     # calculate WFAT data for HOMO (index 0)

# Saving calculated data for future use
MolSaveDataAs!(m, "CO2_molecule.jld2")

# Loading previously saved molecule data
m_ = LoadMolecule("CO2_molecule.jld2")
display(m_)
# [GenericMolecule] Carbon Dioxide
# Asymp coeff of HOMO available
# WFAT data of HOMO available
# #          E (Ha)  occp
# ⋮    ⋮       ⋮      ⋮⋮
# 13 LUMO+1   0.217  ----
# 12 LUMO     0.217  ----
# 11 HOMO    -0.536  -↿⇂-
# 10 HOMO-1  -0.536  -↿⇂-
# 9  HOMO-2  -0.709  -↿⇂-
# 8  HOMO-3  -0.709  -↿⇂-
# ⋮    ⋮       ⋮      ⋮⋮
```

For advanced usages (setting customized grids for quantum chemistry calculations, dealing with non-HOMO orbitals and open-shell molecules, etc.), we suggest reading the [documentation](@ref targets_doc) of the `Targets` module.

-----------------------

### Trajectory Simulation Parameters

The trajectory simulation requires configuring several key parameters.
Below we would guide you to make appropriate choices of the simulation's configuration.
For a comprehensive understanding on the parameters, please refer to the [documentation](@ref traj_simu).

#### Initial Condition and Phase Methods

First of all, the core methods of the trajectory simulation should be selected carefully:
- `init_cond_method`: Which theory to use for generating initial electrons?
- `traj_phase_method`: How quantum effects are handled?

For `init_cond_method`, candidates include:
- SPA-based methods:
  - `ADK`   - Best for small Keldysh parameters (*γ*≪1) when non-adiabatic effects are negligible;
  - `SPANE` - Recommended for most cases, which efficiently includes non-adiabatic effects;
  - `SPA`   - For stronger non-adiabatic scenarios which require full treatment, but much more computational-intensive than `SPANE`.
- `WFAT` - Specialized for molecular targets ([`GenericMolecule`](@ref)) which accounts for the influence of permanent dipoles. However, quantum effects are not available in this method, i.e., only `CTMC` is supported for the phase method.

For `traj_phase_method`:
- `CTMC` - Pure classical trajectories, suitable for scenarios where no inter-cycle interference is involved.
- `QTMC` - Includes basic quantum effects through first-order perturbation.
- `SCTS` - Includes quantum effects beyond perturbation, more accurate than `QTMC`.

#### Detector Configuration

The simulation space and the final momentum grid (i.e., the electron detectors) are controlled by:

- `dimension`: Dimension of the simulation (2 or 3). 2D simulations would be carried out in the ``xy`` plane. For most cases, 2D simulation is sufficient as the electron dynamics is mainly confined in the polarization plane.
- `final_p_max`: Maximum momentum range in each direction. For 2D simulation, specify as a tuple of 2 values `(px_max, py_max)`. For 3D simulation, use a tuple of 3 values `(px_max, py_max, pz_max)`. Values are in atomic units. The momentum ranges from -`p_max` to `p_max`.
- `final_p_num`: Grid resolution in each direction. For 2D, specify as `(nx, ny)` where `nx`, `ny` are number of grid points. For 3D, specify as `(nx, ny, nz)`. Higher resolution provides better details but requires more memory.

Example usage for common scenarios:
```julia
# a 2D example
perform_traj_simulation(
    # ...
    dimension = 2,
    final_p_max = (2.0, 2.0),     # ±2 a.u. in x,y
    final_p_num = (500, 500),     # 500×500 grid
    # ...
)

# a 3D example
perform_traj_simulation(
    # ...
    dimension = 3,
    final_p_max = (2.5, 2.5, 1.0),  # ±2.5 a.u. in x,y, ±1 a.u. in z
    final_p_num = (250, 250, 100),  # 250×250×100 grid
    # ...
)

# The following example performs an effective 2D simulation in x-z plane:
perform_traj_simulation(
    # ...
    dimension = 3,
    final_p_max = (2.5, 0.0, 2.5),
    final_p_num = (400, 1, 400),
    # ...
)
```

!!! tip "Tips"
    A test run with reduced number of trajectories can help determine the optimal `final_p_max`.

#### Sampling Parameters

- Time domain:

```julia
perform_traj_simulation(
    sample_t_intv = (-10.0fs, 10.0fs),  # sampling interval (default in a.u.)
    sample_t_num = 2000,                # number of samples in time domain
    # ...
)
```
Enabling Monte Carlo sampling (`sample_monte_carlo=true`) would generate random samples in the time domain.


- Momentum (``\kkt=(k_\perp,k_z)``) domain:

```julia
# Step sampling:
perform_traj_simulation(
    # ...
    sample_monte_carlo = false,
    ss_kd_max = 2.0,    # transverse momentum range
    ss_kd_num = 200,    # number of kd points
    ss_kz_max = 2.0,    # longitudinal momentum range (3D only)
    ss_kz_num = 200,    # number of kz points (3D only)
    # ...
)

# Monte Carlo sampling:
perform_traj_simulation(
    # ...
    sample_monte_carlo = true,
    mc_kt_num = 100,    # k_t samples per time point
    mc_kd_max = 2.0,    # transverse momentum range
    mc_kz_max = 2.0,    # longitudinal range (3D only)
    # ...
)
```

-----------------------

## Viewing and Analyzing the Results

After running the trajectory simulation, the output file will be generated in the Julia Data Format (JLD2) or HDF5 format. This file contains the photoelectron momentum distribution (PMD) and other necessary information.

To view and analyze the result, you can use the `JLD2.jl` package in Julia. Here is an example of how to open and inspect the output file:

```julia-repl
julia> using JLD2

julia> file = jldopen("ADK-CTMC_4e14_800nm_cos4_2cyc_CP.jld2")
JLDFile .../ADK-CTMC_4e14_800nm_cos4_2cyc_CP.jld2 (read-only)
    ├─🔢 info
    ├─🔢 params_text              # parameters stored in YAML format
    ├─🔢 params                   # parameters stored in a Julia `Dict`
    ├─🔢 px                       # coordinates of the `momentum_spec` on x-axis
    ├─🔢 py                       # coordinates of the `momentum_spec` on y-axis
    ├─🔢 momentum_spec            # PMD data stored in a Julia `Array`
    ├─🔢 ion_prob                 # total ionization probability
    ├─🔢 ion_prob_uncollected     # ionization probability of discarded electrons
    └─🔢 num_effective_traj       # total number of effective trajectories

julia> file["num_effective_traj"]
54240616

julia> file["params_text"] |> print  # equivalent to print(file["params_text"])
init_cond_method: ADK
laser:
  type: Cos4Laser
  peak_int: 4.0e14
  wave_len: 800.0
  cyc_num: 2
  ellip: 1.0
  azi: 0.0
  cep: 0.0
  t_shift: 0.0
target:
  type: HydrogenLikeAtom
  Ip: 0.5
  nucl_charge: 1
  l: 0
  m: 0
  asymp_coeff: 1.0
  soft_core: 1.0e-10
  name: "H"
sample_t_intv: (-100, 100)
sample_t_num: 20000
sample_monte_carlo: false
sample_cutoff_limit: 1.0e-16
traj_phase_method: CTMC
traj_rtol: 1.0e-6
traj_t_final: 120
output_fmt: jld2
output_path: "ADK-CTMC_4e14_800nm_cos4_2cyc_CP.jld2"
final_p_max: (2.5, 2.5)
final_p_num: (500, 500)
ss_kd_max: 2.0
ss_kd_num: 10000
rate_prefix: Full
```

To make plots of the PMD in Julia, you can use the [`Plots.jl`](https://github.com/JuliaPlots/Plots.jl) package.
We included plot scripts in the [`example/`](https://github.com/TheStarAlight/eTraj.jl/tree/master/examples) directory for reference.
Note that running these scripts also requires additional dependencies such as `CodecZlib.jl` and `PyPlot.jl`.
After installing these packages, you can run the plot scripts either in a Julia REPL or directly from the command line, see the ["Running Scripts" section](@ref running-scripts).
If you encounter any issues with plotting, please open an issue in [our GitHub repository](https://github.com/TheStarAlight/eTraj.jl/issues).

!!! note "Guide on Installing the `Plots.jl` & `PyPlot.jl` dependencies"
    The plotting scripts we attached in the `example/` directory requires additional installation of the `Plots.jl` and `PyPlot.jl` packages.
    To install the dependencies, first, follow the instructions [here](@ref config_py) to configure the `PyCall.jl` package and install the `matplotlib` python package.
    Then, install `Plots.jl` and `PyPlot.jl` packages by running the following commands:
    ```julia
    using Pkg
    Pkg.add("Plots")
    Pkg.add("PyPlot")
    ```
