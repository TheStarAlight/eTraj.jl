# <p align="center"> 🎆eTraj.jl </p>

<p align="center">
    Implementation of classical/semiclassical trajectory-based methods in strong-field ionization of atoms and molecules.
</p>

<p align="center">
    <a href="#background">Background</a> •
    <a href="#installation">Installation</a> •
    <a href="#usage">Usage</a> •
    <a href="#examples">Examples</a> •
    <a href="#contributors">Contributors</a> •
    <a href="#license">License</a>
</p>

<!-- <p align="center">
    • Documentation available at [TheStarAlight.github.io/eTraj.jl](https://TheStarAlight.github.io/eTraj.jl) •
</p> -->

<p align="center">
    Last updated: Oct 23, 2024
</p>
<p align="center">
    Version 1.0.0-rc2
</p>

---------------------------

## Background

The interaction between light and matter has attracted widespread interest since the early days of quantum mechanics in the early twentieth century.
With the advent of laser technology, the intensity of light and the precision of spectroscopy has dramatically increased,
which allows us to explore the physics of light-matter interaction under extreme conditions with unprecedented precision and accuracy.
At a high laser intensity above TW/cm², the interaction between light and atoms or molecules can no longer be described by the perturbation theory and a series of novel strong-field phenomena emerges, such as the above-threshold ionization (ATI), tunneling ionization, high-harmonic generation (HHG) and non-sequential double ionization (NSDI).

Theoretical studies of these non-perturbative phenomena have been extensively investigated in the past decades.
Usually, in order to obtain a precise result, a time-dependent Schrödinger equation (TDSE) is solved numerically.
However, solving the TDSE is computationally expensive and resource-demanding, which limits its application to few-dimensional problems.
Moreover, the TDSE is like a black box, and therefore non-transparent for interpretation of the underlying physics.
Apart from TDSE, the strong-field approximation (SFA) is also widely applied to study these problems,
which is based on the assumptions that: (1) the initial state is not affected by the laser field until ionization;
(2) after ionization, the photoelectron is not influenced by the trapping potential (i.e., assuming a short-range potential).
These two approximations simplify the problem to some extent, which allows one to obtain analytical results and unravel the physical pictures of these phenomena.
However, such approximations are not always appropriate for some cases when the Coulomb potential's role becomes significant,
which may lead to incorrect predictions.

To overcome the shortcomings of TDSE, Corkum *et al.* proposed a scheme, where the electron is first ionized from the target through the tunneling mechanism, then acts as a classical electron in the laser field, and finally the photoelectron spectra is obtained by statistics on the ensemble of the electron trajectories [^Corkum_1989] [^Corkum_1993].
This scheme was further developed by Hu *et al.*, in which the initial conditions of the classical electrons and the Coulomb potential of the parent ion are more appropriately taken account [^Hu_1997].
The scheme is named after the Classical-Trajectory Monte-Carlo (CTMC) method, which has been widely adopted.
Further developments of the trajectory-based methods added a phase property to the electron orbits, and utilized the SFA for preparation of the initial conditions,
for example, the Trajectory-based Coulomb-SFA (TC-SFA) [^Yan_2010] [^Yan_2012], the Quantum-Trajectory Monte Carlo (QTMC) [^Li_2014] [^Liu_2016], the Semiclassical Two-Step Model (SCTS) [^Shvetsov-Shilovski_2016] [^Shvetsov-Shilovski_2021],
another approach which is named after the Coulomb Quantum-orbit SFA (CQSFA) [^Lai_2015] [^Maxwell_2017], chose the inverse problem by searching for all trajectories that arrived at the same final momenta.
Such trajectory-based semiclassical methods have a notable advantage over the TDSE and direct SFA methods due to their low demand on computational resources, comparable accuracy to the TDSE, as well as the transparency of the physical picture.

After years of development, various trajectory-based classical/semiclassical methods have flourished, but a unified theoretical framework has yet to be established.
Besides, developing a library that implements existing methods, which is efficient in calculations and is easy to maintain, would greatly facilitate further research on strong-field ionization.
With this aim, we present `eTraj.jl`, a program package written in the Julia language, which provides a general, efficient, and out-of-the-box solution for performing classical/semiclassical trajectory simulations.
Written in a clear and concise manner, this library features versatility, extensibility, and usability.

[^Corkum_1989]: P. B. Corkum, N. H. Burnett, and F. Brunel, Above-threshold ionization in the long-wavelength limit, *Phys. Rev. Lett.* **62**, 1259 (1989). DOI: [10.1103/PhysRevLett.62.1259](https://dx.doi.org/10.1103/PhysRevLett.62.1259)

[^Corkum_1993]: P. B. Corkum, Plasma perspective on strong field multiphoton ionization, *Phys. Rev. Lett.* **71**, 1994 (1993). DOI: [10.1103/PhysRevLett.71.1994](https://doi.org/10.1103/PhysRevLett.71.1994)

[^Hu_1997]: B. Hu, J. Liu, and S.-g. Chen, Plateau in Above-Threshold-Ionization Spectra and Chaotic Behavior in Rescattering Processes. *Phys. Lett. A* **236**, 533–542 (1997). DOI: [10.1016/S0375-9601(97)00811-6](https://dx.doi.org/10.1016/S0375-9601(97)00811-6)

[^Yan_2010]: T.-M. Yan, S. V. Popruzhenko, M. J. J. Vrakking, and D. Bauer, Low-energy structures in strong field ionization revealed by quantum orbits, *Phys. Rev. Lett.* **105**, 253002 (2010). DOI: [10.1103/PhysRevLett.105.253002](https://doi.org/10.1103/PhysRevLett.105.253002)

[^Yan_2012]: T.-M. Yan and D. Bauer, Sub-barrier coulomb effects on the interference pattern in tunneling-ionization photoelectron spectra, *Phys. Rev. A* **86**, 53403 (2012). DOI: [10.1103/PhysRevA.86.053403](https://doi.org/10.1103/PhysRevA.86.053403)

[^Li_2014]: M. Li, J.-W. Geng, H. Liu, Y. Deng, C. Wu, L.-Y. Peng, Q. Gong, and Y. Liu, Classical-quantum correspondence for above-threshold ionization, *Phys. Rev. Lett.* **112**, 113002 (2014). DOI: [10.1103/PhysRevLett.112.113002](https://doi.org/10.1103/PhysRevLett.112.113002)

[^Liu_2016]: M.-M. Liu, M. Li, C. Wu, Q. Gong, A. Staudte, and Y. Liu, Phase structure of strong-field tunneling wave packets from molecules, *Phys. Rev. Lett.* **116**, 163004 (2016). DOI: [10.1103/PhysRevLett.116.163004](https://doi.org/10.1103/PhysRevLett.116.163004)

[^Shvetsov-Shilovski_2016]: N. I. Shvetsov-Shilovski, M. Lein, L. B. Madsen, E. Räsänen, C. Lemell, J. Burgdörfer, D. G. Arbó, and K. Tőkési, Semiclassical two-step model for strong-field ionization, *Phys. Rev. A* **94**, 013415 (2016). DOI: [10.1103/PhysRevA.94.013415](https://doi.org/10.1103/PhysRevA.94.013415)

[^Shvetsov-Shilovski_2021]: N. I. Shvetsov-Shilovski, Semiclassical two-step model for ionization by a strong laser pulse: Further developments and applications, *Eur. Phys. J. D* **75**, 130 (2021). DOI: [10.1140/epjd/s10053-021-00095-8](https://doi.org/10.1140/epjd/s10053-021-00095-8)

[^Lai_2015]: X.-Y. Lai, C. Poli, H. Schomerus, and C. F. D. M. Faria, Influence of the coulomb potential on above-threshold ionization: A quantum-orbit analysis beyond the strong-field approximation, *Phys. Rev. A* **92**, 043407 (2015). DOI: [10.1103/PhysRevA.92.043407](https://doi.org/10.1103/PhysRevA.92.043407)

[^Maxwell_2017]: A. S. Maxwell, A. Al-Jawahiry, T. Das, and C. F. D. M. Faria, Coulomb-corrected quantum interference in above-threshold ionization: Working towards multitrajectory electron holography, *Phys. Rev. A* **96**, 023420 (2017). DOI: [10.1103/PhysRevA.96.023420](https://doi.org/10.1103/PhysRevA.96.023420)

---------------------------

## Installation

### Prerequisites

- *Minimum prerequisites* : Julia ≥ 1.9

- *MO-ADK/MO-SFA and WFAT molecular calculations* :
    Data for some small molecules are available in the molecule database.
    If the user wants to perform molecular calculations with customized parameters, the platform should be **Linux or macOS**, and having *Python 3* with the [PySCF](https://github.com/pyscf/pyscf) python package installed and the [PyCall.jl](https://github.com/JuliaPy/PyCall.jl) julia package successfully built.

### Installing the package

This package is currently not in julia's general registry, but can be added through the repository URL:

```julia
using Pkg
Pkg.add("https://github.com/TheStarAlight/eTraj.jl.git")
# In pkg mode of REPL:
# (@v1.9) pkg> add https://github.com/TheStarAlight/eTraj.jl.git
```

To enter the pkg mode of REPL, type `]` in REPL, and the `pkg>` prompt will appear, replacing the `julia>`.

For offline installation:

```julia
using Pkg
Pkg.add(url="/path/to/eTraj/")
# In pkg mode of REPL:
# (@v1.9) pkg> add /path/to/eTraj/
```

It is suggested to test the package to check if the functions (especially molecular calculations) work on your platform:

```julia
Pkg.test("eTraj")
# In pkg mode of REPL:
# (@v1.9) pkg> test eTraj
```

### Configuring Python and PySCF

Currently, the calculation of molecules' asymptotic coefficients (for MO-ADK/MO-SFA) and WFAT coefficients rely on the [`PySCF`](https://github.com/pyscf/pyscf) python package. eTraj calls the PySCF using the [`PyCall.jl`](https://github.com/JuliaPy/PyCall.jl) package.

There are two ways to set up the Python environment used by `PyCall`:

1. using your local Python environment by specifying the path of your Python executable in `ENV["PYTHON"]` and build the PyCall package.
2. using a private Python environment managed by the [`Conda.jl`](https://github.com/JuliaPy/Conda.jl), which is implicitly installed by the `PyCall` package by default;

#### Using the local Python environment

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
pip install pyscf
```

#### Using Conda.jl

Before installing `eTraj`, install `Conda` first:

```julia
using Pkg
Pkg.add("Conda")
```

Then call `pip` within `Conda` to install `PySCF`:

```julia
using Conda
Conda.pip_interop(true)
Conda.pip("install", "pyscf")
```

**Note**:
Since the `PySCF` does not support the Windows platform, the molecular calculation must be performed on a Linux or macOS platform.
However, for Windows users, they may install the [WSL](https://learn.microsoft.com/en-us/windows/wsl/install) (Windows Subsystem for Linux), which supports the `PySCF`.

## Usage

The usage of the program is simple and straightforward.
The `eTraj.perform_traj_simulation` method serves as a public entrance to performing a trajectory simulation.
The method would automatically detect number of available threads (specified by passing command-line arguments `-t <thread_num>` when starting julia) and run the trajectory simulation in parallel.

### Parameters of the `perform_traj_simulation` method

#### Required parameters

- `init_cond_method`    : Method used to determine the initial conditions of electrons.
  - Candidates: `:ADK`, `:SPA` (SFA-SPA), `:SPANE` (SFA-SPANE) for targets of type `SAEAtomBase` or `MoleculeBase`; `:WFAT` for `MoleculeBase` targets.
- `laser::Laser`        : Parameters of the laser field. See the [`Laser`](#the-lasers-module) module for details.
- `target::Target`      : Parameters of the target. See the [`Targets`](#the-targets-module) module for details.
- `dimension = 2|3`     : Dimensionality of simulation.
  - 2D simulation is carried out in the xy plane.
- `sample_t_intv`       : Time interval for sampling initial electrons.
  - Format: `(start,stop)`
  - Unit: pass numerically in **a.u.** or pass as a `Unitful.Quantity`.
- `sample_t_num`        : Number of time samples.
- `traj_t_final`        : Final time of each trajectory simulation
  - Unit: numerically in **a.u.** or pass as a `Unitful.Quantity`.
- `final_p_max`         : Boundaries of final momentum grid. Grid ranges from `-pxMax` to `+pxMax` in the x direction, and the same for y and z directions.
  - Format: `(pxMax,pyMax[,pzMax])`
- `final_p_num`         : Numbers of final momentum grid points. If a value is `1`, electrons will be collected regardless of the momentum on that dimension.
  - Format: `(pxNum,pyNum[,pzNum])`

#### Required parameters for step-sampling methods

- `ss_kd_max`   : Boundary of k⟂ samples (in a.u.). k⟂ ranges from `-ss_kd_max` to `+ss_kd_max`.
- `ss_kd_num`   : Number of k⟂ samples (in a.u.).
- `ss_kz_max`   : [3D only] Boundary of kz samples (in a.u.). kz ranges from `-ss_kz_max` to `+ss_kz_max`.
- `ss_kz_num`   : [3D only] Number of kz samples (in a.u.).

#### Required parameters for Monte-Carlo sampling methods

- `mc_kt_num`   : Number of kt samples in a single time sample.
- `mc_kd_max`   : Boundary of kd. kd ranges from `-mc_kd_max` to `+mc_kd_max`.
- `mc_kz_max`   : [3D only] Boundary of kz. kz ranges from `-mc_kz_max` to `+mc_kz_max`.

#### Optional parameters

- `traj_phase_method`   : Method used to determine classical trajectories' phase.
  - Candidates: `:CTMC` (default), `:QTMC`, and `:SCTS`.
- `traj_rtol`           : Relative error tolerance for solving classical trajectories (default `1e-6`).
- `output_fmt`          : Output file format.
  - Candidates: `:jld2` (JLD2, default) and `:h5` (HDF5).
- `output_compress`     : Determines whether output files are compressed or not (default `true`).
  - Note: For JLD2 output format, compression requires explicit installation of the `CodecZlib` package.
- `output_path`         : Path to output file.
- `sample_cutoff_limit` : Probability cutoff limit for sampled electrons (default `1e-16`). Electrons with probabilities lower than the limit would be discarded.
- `sample_monte_carlo`  : Determines whether Monte-Carlo sampling is used when generating electron samples (default `false`).

#### Optional parameter for atomic SFA-SPA, SFA-SPANE and ADK methods

- `rate_prefix` : Prefix of the exponential term in the ionization rate (default `:Full`).
  - `:Exp` indicates no prefix; `:Pre` and `:PreCC` indicates inclusion of the prefactor with/without Coulomb correction; `:Jac` indicates inclusion of the Jacobian factor which is related to the sampling method; `:Full` is equivalent to `Set([:PreCC,:Jac])`.
      To combine `:Pre` and `:Jac`, pass `Set([:Pre,:Jac])`.

#### Optional parameter for target `MoleculeBase`

- `mol_orbit_ridx`  : Index of selected orbital relative to the HOMO (e.g., `0` indicates HOMO, and `-1` indicates HOMO-1.)
  - For open-shell molecules, according to α/β spins, should be passed in format `(spin, idx)` where for α orbitals spin=`1` and for β orbitals spin=`2`.

#### Optional parameter for interface

- `show_progress`   : Whether to display progress bar (default `true`).

### The `Lasers` module

The `Lasers` module provides abstraction of laser pulses, including the following types:

```text
⊚ Laser
├ ⊚ MonochromaticLaser
│ ├ ⦶ Cos2Laser
│ ├ ⦶ Cos4Laser
│ └ ⦶ GaussianLaser
└ ⦶ BichromaticLaser
```

where `⊚` indicates an abstract type, `⦶` indicates a concrete type and the tree structure shows the inheritance relation.

The initializer of each type is defined as follows:

```julia
Cos4Laser(peak_int, wave_len|ang_freq, cyc_num|duration, ellip [,azi=0] [,cep=0] [,t_shift=0])
Cos2Laser(peak_int, wave_len|ang_freq, cyc_num|duration, ellip [,azi=0] [,cep=0] [,t_shift=0])
GaussianLaser(peak_int, wave_len|ang_freq, spread_cyc_num|spread_duration|FWHM_duration, ellip [,azi=0] [,cep=0] [,t_shift=0])
BichromaticLaser(l1::MonochromaticLaser, l2::MonochromaticLaser [,delay=0])
```

for the detailed description of each parameter, please refer to [the documentation](https://TheStarAlight.github.io/eTraj.jl/).

Some examples:

```julia
julia> using eTraj.Lasers

# without units, following default units
julia> l = Cos2Laser(peak_int=4e14, wave_len=800.0, cyc_num=2.0, ellip=0.0)
[MonochromaticLaser] Envelope cos², peak intensity 4.0e+14 W/cm², wavelen=800 nm, 2 cycle(s), ε=0 [linearly polarized]

# import the Units submodule to access units
julia> using eTraj.Units

julia> l = Cos4Laser(peak_int=0.4PW/cm^2, ang_freq=1.5498eV, duration=26.7fs, ellip=1.0, cep=90°)
[MonochromaticLaser] Envelope cos⁴, peak intensity 4.0e+14 W/cm², wavelen=800.00 nm, 10.01 cycle(s), ε=1 [circularly polarized], CEP=0.50 π

julia> l = GaussianLaser(peak_int=4e14W/cm^2, wave_len=400nm, FWHM_duration=20fs, ellip=-1)
[MonochromaticLaser] Envelope Gaussian, peak intensity 4.0e+14 W/cm², wavelen=400 nm, temporal width 9.00 cycle(s) [FWHM 20.00 fs], ε=-1 [circularly polarized]

julia> l = BichromaticLaser(l1=Cos4Laser(peak_int=1.0PW/cm^2, wave_len=800nm, cyc_num=10, ellip=1), l2=Cos4Laser(peak_int=1.0PW/cm^2, wave_len=400nm, cyc_num=20, ellip=-1), delay=0.5fs)
[BichromaticLaser] delay Δt = 20.67 a.u. (0.50 fs)
├ [MonochromaticLaser] Envelope cos⁴, peak intensity 1.0e+15 W/cm², wavelen=800 nm, 10 cycle(s), ε=1 [circularly polarized]
└ [MonochromaticLaser] Envelope cos⁴, peak intensity 1.0e+15 W/cm², wavelen=400 nm, 20 cycle(s), ε=-1 [circularly polarized]
```

### The `Targets` module

The `Targets` module provides abstraction of targets, including the following types:

```text
⊚ Target
├ ⊚ SAEAtomBase
│ ├ ⦶ HydrogenLikeAtom
│ └ ⦶ SAEAtom
└ ⊚ MoleculeBase
  └ ⦶ GenericMolecule
```

The initializer of each type is defined as follows:

```julia
HydrogenLikeAtom(Ip, Z [,l=0] [,m=0] [,asymp_coeff=:hartree|<coeff>] [,quan_ax_θ=0] [,quan_ax_ϕ=0] [,soft_core=1e-10] [,name])
SAEAtom(Ip, Z [,l=0] [,m=0] [,asymp_coeff=:hartree|<coeff>] [,quan_ax_θ=0] [,quan_ax_ϕ=0] [,a1,b1,a2,b2,a3,b3] [,soft_core=1e-10] [,name])
GenericMolecule(atoms, atom_coords [,charge=0] [,spin=0] [,name] [,rot_α=0] [,rot_β=0] [,rot_γ=0])
# this method loads a GenericMolecule from an external file.
LoadMolecule(ext_data_path; [rot_α=0] [,rot_β=0] [,rot_γ=0])
```

Some examples of using the `HydrogenLikeAtom` and `SAEAtom` types:

```julia
julia> using eTraj.Targets

julia> t = HydrogenLikeAtom(Ip=0.5, Z=1, name="H")
[HydrogenLikeAtom] Atom H, Ip=0.5000 (13.61 eV), Z=1

julia> using eTraj.Units

julia> t = SAEAtom(Ip=12.13eV, Z=1, l=1, a1=51.35554, b1=2.111554, a2=-99.92747, b2=3.737221, a3=1.644457, b3=0.4306465, asymp_coeff=1.3, name="Xe")
[SAEAtom] Atom Xe (p orbital, m=0), Ip=0.4458 (12.13 eV), Z=1
```

Some examples of using the `GenericMolecule` type:

```julia
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

We also provide some preset atoms and molecules for instant use:

```julia
julia> using eTraj.Targets, eTraj.Units

# example of get_atom
julia> t = get_atom("He1p")
[HydrogenLikeAtom] Atom He⁺, Ip=1.0000 (27.21 eV), Z=2

julia> t = get_atom("Xe"; m=1, quan_ax_θ=90°, quan_ax_ϕ=0°) # other parameters can be passed via keyword arguments
[SAEAtom] Atom Xe (p orbital, m=1), Ip=0.4458 (12.13 eV), Z=1, θϕ=(90.0°,0.0°)

# example of get_mol
julia> mol = get_mol("Nitric Oxide")
[GenericMolecule] Nitric Oxide (NO)
Asymp coeff of α-HOMO & α-LUMO available
WFAT data of α-HOMO & α-LUMO available
#          Eα(Ha)  occp  Eβ(Ha)
⋮    ⋮        ⋮     ⋮⋮      ⋮     ⋮
10 LUMO+1   0.378  ----   0.386 LUMO+2
9  LUMO    -0.415  ----   0.156 LUMO+1
8  HOMO    -0.415  -↿--   0.156 LUMO
7  HOMO-1  -0.696  -↿⇂-  -0.613 HOMO
6  HOMO-2  -0.784  -↿⇂-  -0.613 HOMO-1
5  HOMO-3  -0.784  -↿⇂-  -0.654 HOMO-2
4  HOMO-4  -0.969  -↿⇂-  -0.884 HOMO-3
⋮    ⋮        ⋮     ⋮⋮      ⋮     ⋮
```

The available keys can be accessed by invoking `get_available_atoms()` and `get_available_mols()`.

---------------------------

## Examples

The corresponding simulation and plotting scripts can be found in the [`examples/`](https://github.com/TheStarAlight/eTraj.jl/tree/master/examples) directory.

<!-- For detailed descriptions, see the [documentation](). -->

### Attoclock Experiment using different initial condition methods

```julia
using eTraj
using eTraj.Targets, eTraj.Lasers, eTraj.Units

l = Cos4Laser(peak_int=0.4PW/cm^2, wave_len=800.0nm, cyc_num=2, ellip=1.0)
t = get_atom("H")

for init_cond in [:ADK, :SPANE, :SPA]
    perform_traj_simulation(
        init_cond_method    = init_cond,
        laser               = l,
        target              = t,
        dimension           = 2,            # 2D simulation, x-y plane only
        sample_t_intv       = (-100,100),   # equivalent to `(-2.42fs, 2.42fs)`
        sample_t_num        = 20000,        # will sample 20000 equidistant time points between -100 and 100 a.u.
        traj_t_final        = 120,          # the traj end at 120 a.u., equivalent to `2.90fs`
        final_p_max         = (2.5,2.5),    # the momentum spec collection grid's border (-2.5 to +2.5 a.u.)
        final_p_num         = (500,500),    # the momentum spec collection grid's size (500x500)
        ss_kd_max           = 2.0,
        ss_kd_num           = 10000,        # will sample 10000 equidistant k⟂ points between -2 to +2 a.u.
        output_path         = "$(init_cond)-CTMC_4e14_800nm_cos4_2cyc_CP.jld2",
        traj_phase_method   = :CTMC
    )
end
```

### Interaction with linearly-polarized pulses and the influence of phase methods

```julia
# 8-cycle
using eTraj
using eTraj.Targets, eTraj.Lasers, eTraj.Units

l = Cos2Laser(peak_int=90.0TW/cm^2, wave_len=800.0nm, cyc_num=8, ellip=0.0)
t = get_atom("H")

for phase_method in [:QTMC, :SCTS]
    perform_traj_simulation(
        init_cond_method    = :ADK,
        laser               = l,
        target              = t,
        dimension           = 3,
        sample_t_intv       = (-350,350),
        sample_t_num        = 50000,
        traj_t_final        = 500,
        # == perform the simulation only in the x-z plane ==
        final_p_max         = (1.0,0.0,1.0),
        final_p_num         = (500,1,500),
        ss_kd_max           = 0.0,
        ss_kd_num           = 1,
        # ==================================================
        ss_kz_max           = 1.0,
        ss_kz_num           = 20000,
        output_path         = "ADK-$(phase_method)_9e13_800nm_8cyc_LP_ExpRate.jld2",
        traj_phase_method   = phase_method,
        rate_prefix         = :Exp
    )
end
```

```julia
# 1-cycle
using eTraj
using eTraj.Targets, eTraj.Lasers, eTraj.Units

l = Cos2Laser(peak_int=90.0TW/cm^2, wave_len=800.0nm, cyc_num=1, cep=π/2, ellip=0.0)
t = get_atom("H"; soft_core=1e-12)

for phase_method in [:QTMC, :SCTS]
    perform_traj_simulation(
        init_cond_method    = :ADK,
        laser               = l,
        target              = t,
        dimension           = 3,
        sample_t_intv       = (-50,50),
        sample_t_num        = 30000,
        traj_t_final        = 100,
        # == perform the simulation only in the x-z plane ==
        final_p_max         = (1.0,0.0,1.0),
        final_p_num         = (500,1,500),
        ss_kd_max           = 0,
        ss_kd_num           = 1,
        # ==================================================
        ss_kz_max           = 1.5,
        ss_kz_num           = 10000,
        output_path         = "ADK-$(phase_method)_9e13_800nm_1cyc_LP_ExpRate.jld2",
        traj_phase_method   = phase_method,
        rate_prefix         = :Exp
    )
end
```

### Interaction with an $\omega-2\omega$ bichromatic clover-shaped laser

```julia
using eTraj
using eTraj.Targets, eTraj.Lasers, eTraj.Units

for int in [1e14, 3e14, 5e14, 7e14]
    @info "Running I0=$(int) W/cm^2"
    l1 = Cos2Laser(peak_int=int*W/cm^2, wave_len=800.0nm, cyc_num=8,  ellip= 1.0)
    l2 = Cos2Laser(peak_int=int*W/cm^2, wave_len=400.0nm, cyc_num=16, ellip=-1.0)
    l = BichromaticLaser(l1=l1, l2=l2)
    t = get_atom("H")
    perform_traj_simulation(
        init_cond_method    = :ADK,
        laser               = l,
        target              = t,
        dimension           = 2,
        sample_t_intv       = (-350,350),
        sample_t_num        = 10000,
        traj_t_final        = 450,
        final_p_max         = (2.5,2.5),
        final_p_num         = (500,500),
        ss_kd_max           = 1.0,
        ss_kd_num           = 5000,
        output_path         = "ADK-SCTS_Bichromatic_$(int)_800+400nm_8+16cycs_CounterCP.jld2",
        traj_phase_method   = :SCTS,
        rate_prefix         = Set([:Pre,:Jac])  # Coulomb correction unavailable for bichromatic laser
    )
end
```

### WFAT-CTMC simulation of molecular targets

```julia
using eTraj
using eTraj.Targets, eTraj.Lasers, eTraj.Units

l = Cos2Laser(peak_int=4e14W/cm^2, wave_len=800.0nm, cyc_num=6, ellip=1.0)
t = [get_mol("Hydrogen"; rot_β=90°),
     get_mol("Carbon Monoxide"; rot_β=90°),
     get_mol("Oxygen"; rot_β=90°),
     get_mol("Oxygen"; rot_β=90°),
     get_mol("Benzene"; rot_β=90°),
     get_mol("Benzene"; rot_β=90°)]
orbit_ridx = [0, 0, (1,0), (1,-1), 0, -1]
path = [
    "WFAT-CTMC_Hydrogen_HOMO_4e14_800nm_6cyc_CP.jld2",
    "WFAT-CTMC_CarbonMonoxide_HOMO_4e14_800nm_6cyc_CP.jld2",
    "WFAT-CTMC_Oxygen_α-HOMO_4e14_800nm_6cyc_CP.jld2",
    "WFAT-CTMC_Oxygen_α-HOMO-1_4e14_800nm_6cyc_CP.jld2",
    "WFAT-CTMC_Benzene_HOMO_4e14_800nm_6cyc_CP.jld2",
    "WFAT-CTMC_Benzene_HOMO-1_4e14_800nm_6cyc_CP.jld2"
]
for i in eachindex(t)
    perform_traj_simulation(
        init_cond_method    = :WFAT,
        laser               = l,
        target              = t[i],
        dimension           = 2,
        sample_t_intv       = (-300,300),
        sample_t_num        = 10000,
        traj_t_final        = 350,
        final_p_max         = (2.0,2.0),
        final_p_num         = (500,500),
        ss_kd_max           = 2.0,
        ss_kd_num           = 5000,
        output_path         = path[i],
        traj_phase_method   = :CTMC,    # WFAT supports CTMC only
        mol_orbit_ridx      = orbit_ridx[i]
    )
end
```

---------------------------

## Contributors

- [Mingyu Zhu](https://github.com/TheStarAlight) @ ECNU
- [Hongcheng Ni](https://faculty.ecnu.edu.cn/_s29/nhc_en/main.psp) @ ECNU

---------------------------

## License

This package is licensed under the Apache 2.0 license, and is copyrighted by Mingyu Zhu and the other contributors.
