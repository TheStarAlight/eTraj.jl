# 🎆eTraj.jl

*Implementation of classical/semiclassical trajectory-based methods in strong-field ionization of atoms and molecules.*

---------------------------

```@contents
Pages = ["index.md"]
```

---------------------------

## Background

The interaction between light and matter has attracted widespread interest since the early days of quantum mechanics.
With the advent of laser technology, the intensity of light and the precision of spectroscopy has dramatically increased,
which allows us to explore the physics of light-matter interaction under extreme conditions with unprecedented precision and accuracy.
At a high laser intensity above TW/cm², the interaction between light and atoms or molecules can no longer be described by the perturbation theory and a series of novel strong-field phenomena emerges, such as the above-threshold ionization (ATI), tunneling ionization, high-harmonic generation (HHG) and non-sequential double ionization (NSDI).



Theoretical studies of these non-perturbative phenomena have been extensively investigated in the past decades.
Usually, in order to obtain a precise result, a time-dependent Schrödinger equation (TDSE) is solved numerically.
However, solving the TDSE is computationally expensive and resource-demanding, which limits its application to few-dimensional problems.
Moreover, the TDSE is like a black box, offering limited transparency for interpreting the underlying physics.
Apart from TDSE, the strong-field approximation (SFA) is also widely applied to study these problems,
which is based on the assumptions that: (1) the initial state is not affected by the laser field until ionization;
(2) after ionization, the photoelectron is not influenced by the trapping potential (i.e., assuming a short-range potential).
These two approximations simplify the problem, which allows one to obtain analytical results and unravel the physical pictures of these phenomena.
However, such approximations are not always applicable, especially when the Coulomb potential's role becomes significant,
which may lead to incorrect predictions.

To address these limitations, the scheme of *Classical-Trajectory Monte-Carlo* (CTMC) method [^Abrines_1966] [^Olson_1977] can be adopted, where a microcanonical ensemble of classical electrons is prepared and evolved under the laser interaction or charged-particle impact.
This scheme has been further developed to account for the initial stage of tunneling ionization by setting the initial conditions of the classical electrons at the tunnel exit [^Corkum_1989] [^Corkum_1993] [^Hu_1997].
The final photoelectron momentum distribution (PMD) is obtained through statistical analysis of the electron trajectories.
While CTMC relies on purely classical electron trajectories, quantum effects can be largely retained by incorporating a phase into the electron trajectories.
Examples include  the Trajectory-based Coulomb-SFA (TC-SFA) [^Yan_2010] [^Yan_2012], the Quantum-Trajectory Monte Carlo (QTMC) [^Li_2014] [^Liu_2016], and the Semiclassical Two-Step Model (SCTS) [^ShvetsovShilovski_2016] [^ShvetsovShilovski_2021].
Another approach, the Coulomb Quantum-orbit SFA (CQSFA) [^Lai_2015] [^Maxwell_2017], addresses the inverse problem by identifying all trajectories that result in the same final momenta.
These trajectory-based semiclassical methods offer notable advantages over the TDSE and direct SFA methods due to their lower demand on computational resources, as well as the clarity they provide in understanding the physical picture.

After years of development, various trajectory-based classical/semiclassical methods have emerged, yet a unified theoretical framework remains to be established.
Besides, developing a library that implements existing methods, which is efficient in calculations and is easy to maintain, would greatly facilitate further research on strong-field ionization.
With this goal in mind, we introduce `eTraj.jl`, a program package written in the Julia language, which provides a general, efficient, and out-of-the-box solution for performing classical/semiclassical trajectory simulations.
This library is written in a clear and concise manner, ensuring versatility, extensibility, and usability.

[^Abrines_1966]: R. Abrines and I. C. Percival, Classical theory of charge transfer and ionization of hydrogen atoms by protons, *Proc. Phys. Soc.* **88**, 861 (1966). DOI: [10.1088/0370-1328/88/4/306](https://doi.org/10.1088/0370-1328/88/4/306)

[^Olson_1977]: R. E. Olson and A. Salop, Charge-transfer and impact-ionization cross sections for fully and partially stripped positive ions colliding with atomic hydrogen, *Phys. Rev. A* **16**, 531 (1977). DOI: [10.1103/PhysRevA.16.531](https://doi.org/10.1103/PhysRevA.16.531)

[^Corkum_1989]: P. B. Corkum, N. H. Burnett, and F. Brunel, Above-threshold ionization in the long-wavelength limit, *Phys. Rev. Lett.* **62**, 1259 (1989). DOI: [10.1103/PhysRevLett.62.1259](https://dx.doi.org/10.1103/PhysRevLett.62.1259)

[^Corkum_1993]: P. B. Corkum, Plasma perspective on strong field multiphoton ionization, *Phys. Rev. Lett.* **71**, 1994 (1993). DOI: [10.1103/PhysRevLett.71.1994](https://doi.org/10.1103/PhysRevLett.71.1994)

[^Hu_1997]: B. Hu, J. Liu, and S.-g. Chen, Plateau in Above-Threshold-Ionization Spectra and Chaotic Behavior in Rescattering Processes. *Phys. Lett. A* **236**, 533–542 (1997). DOI: [10.1016/S0375-9601(97)00811-6](https://dx.doi.org/10.1016/S0375-9601(97)00811-6)

[^Yan_2010]: T.-M. Yan, S. V. Popruzhenko, M. J. J. Vrakking, and D. Bauer, Low-energy structures in strong field ionization revealed by quantum orbits, *Phys. Rev. Lett.* **105**, 253002 (2010). DOI: [10.1103/PhysRevLett.105.253002](https://doi.org/10.1103/PhysRevLett.105.253002)

[^Yan_2012]: T.-M. Yan and D. Bauer, Sub-barrier coulomb effects on the interference pattern in tunneling-ionization photoelectron spectra, *Phys. Rev. A* **86**, 53403 (2012). DOI: [10.1103/PhysRevA.86.053403](https://doi.org/10.1103/PhysRevA.86.053403)

[^Li_2014]: M. Li, J.-W. Geng, H. Liu, Y. Deng, C. Wu, L.-Y. Peng, Q. Gong, and Y. Liu, Classical-quantum correspondence for above-threshold ionization, *Phys. Rev. Lett.* **112**, 113002 (2014). DOI: [10.1103/PhysRevLett.112.113002](https://doi.org/10.1103/PhysRevLett.112.113002)

[^Liu_2016]: M.-M. Liu, M. Li, C. Wu, Q. Gong, A. Staudte, and Y. Liu, Phase structure of strong-field tunneling wave packets from molecules, *Phys. Rev. Lett.* **116**, 163004 (2016). DOI: [10.1103/PhysRevLett.116.163004](https://doi.org/10.1103/PhysRevLett.116.163004)

[^ShvetsovShilovski_2016]: N. I. Shvetsov-Shilovski, M. Lein, L. B. Madsen, E. Räsänen, C. Lemell, J. Burgdörfer, D. G. Arbó, and K. Tőkési, Semiclassical two-step model for strong-field ionization, *Phys. Rev. A* **94**, 013415 (2016). DOI: [10.1103/PhysRevA.94.013415](https://doi.org/10.1103/PhysRevA.94.013415)

[^ShvetsovShilovski_2021]: N. I. Shvetsov-Shilovski, Semiclassical two-step model for ionization by a strong laser pulse: Further developments and applications, *Eur. Phys. J. D* **75**, 130 (2021). DOI: [10.1140/epjd/s10053-021-00095-8](https://doi.org/10.1140/epjd/s10053-021-00095-8)

[^Lai_2015]: X.-Y. Lai, C. Poli, H. Schomerus, and C. F. D. M. Faria, Influence of the coulomb potential on above-threshold ionization: A quantum-orbit analysis beyond the strong-field approximation, *Phys. Rev. A* **92**, 043407 (2015). DOI: [10.1103/PhysRevA.92.043407](https://doi.org/10.1103/PhysRevA.92.043407)

[^Maxwell_2017]: A. S. Maxwell, A. Al-Jawahiry, T. Das, and C. F. D. M. Faria, Coulomb-corrected quantum interference in above-threshold ionization: Working towards multitrajectory electron holography, *Phys. Rev. A* **96**, 023420 (2017). DOI: [10.1103/PhysRevA.96.023420](https://doi.org/10.1103/PhysRevA.96.023420)

---------------------------

## Installation

### Prerequisites

- *Minimum prerequisites* : Julia ≥ 1.9

- *MO-ADK/MO-SFA and WFAT molecular calculations* :
    Data for some small molecules are available in the molecule database.
    If the user wants to perform molecular calculations with customized parameters, the platform should be **Linux or macOS**, and having *Python 3* with the [`PySCF`](https://github.com/pyscf/pyscf) python package installed and the [`PyCall.jl`](https://github.com/JuliaPy/PyCall.jl) julia package successfully built.

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
pip install pyscf==2.3.0
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
Conda.pip("install", "pyscf==2.3.0")
```

**Note**:
Since the `PySCF` does not support the Windows platform, the molecular calculation must be performed on a Linux or macOS platform.
However, for Windows users, they may install the [WSL](https://learn.microsoft.com/en-us/windows/wsl/install) (Windows Subsystem for Linux), which supports the `PySCF`.

---------------------------

## Contributors

- [Mingyu Zhu](https://github.com/TheStarAlight) @ ECNU
- [Hongcheng Ni](https://faculty.ecnu.edu.cn/_s29/nhc_en/main.psp) @ ECNU

---------------------------

## License

This package is licensed under the Apache 2.0 license, and is copyrighted by Mingyu Zhu and the other contributors.
