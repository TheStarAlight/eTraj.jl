# Example: Overview

This section presents selected examples demonstrating key applications of `eTraj`, illustrating its capabilities across representative scenarios.
The example and the corresponding post-processing codes are provided in the `examples/` directory.

The implementation consists of two types of scripts: simulation scripts (named `test_*.jl`) that execute trajectory simulations and generate JLD2 output files, and visualization scripts (named `plot_*.jl`) that process these files to produce figures.
To run the plotting scripts, the user must install and successfully build the `Plots.jl` and `PyPlot.jl` packages.

The execution time of the program is primarily determined by the number of *effective* trajectories, defined as those not filtered out due to probabilities below the `sample_cutoff_limit` threshold. Additional factors influencing computational efficiency include the choice of initial condition and phase methods, the system's dimensionality (specified by the `dimension` parameter), the relative error tolerance of the ODE integrator (specified by the `traj_rtol` parameter), and the occurrence of electron recollisions, which necessitate finer time steps for accurate resolution.
A summary of the execution times for the presented examples is provided in the table below.

The execution time does not include the package loading or precompilation time.
The examples are run on an AMD Ryzen 9 7950X CPU with 12 available cores and 32 GB of RAM, on WSL Ubuntu 22.04 LTS, Julia version 1.11.1.

| Example                   | Configuration      | Cores | Eff. traj.  | Exec. time | Speed        |
|---------------------------|--------------------|-------|-------------|------------|--------------|
| `test_2cycs_CP`           | ADK-CTMC           | 12    | 54.2M       | 3'53''     | 4.3 μs/traj  |
|                           | ADK-CTMC           | 8     | 54.2M       | 4'43''     | 5.2 μs/traj  |
|                           | ADK-CTMC           | 6     | 54.2M       | 5'32''     | 6.1 μs/traj  |
|                           | ADK-CTMC           | 4     | 54.2M       | 7'23''     | 8.2 μs/traj  |
|                           | ADK-CTMC           | 2     | 54.2M       | 11'22''    | 12.5 μs/traj |
|                           | ADK-CTMC           | 1     | 54.2M       | 20'49''    | 23.0 μs/traj |
|                           | SPANE-CTMC         | 12    | 69.7M       | 4'25''     | 3.8 μs/traj  |
|                           | SPA-CTMC           | 12    | 83.9M       | 22'38''    | 16.2 μs/traj |
| `test_8cycs_CP`           | ADK-CTMC           | 12    | 131M        | 15'20''    | 7.0 μs/traj  |
|                           | ADK-QTMC           | 12    | 131M        | 17'59''    | 8.2 μs/traj  |
|                           | ADK-SCTS           | 12    | 131M        | 20'32''    | 9.4 μs/traj  |
| `test_1cyc_CP`            | ADK-CTMC           | 12    | 22.8M       | 1'42''     | 4.5 μs/traj  |
|                           | ADK-QTMC           | 12    | 22.8M       | 2'00''     | 5.3 μs/traj  |
|                           | ADK-SCTS           | 12    | 22.8M       | 2'02''     | 5.4 μs/traj  |
| `test_Bichromatic_CCP`    | I₀ = 1×10¹⁴ W/cm²  | 12    | 18.7M       | 3'54''     | 12.5 μs/traj |
|                           | I₀ = 3×10¹⁴ W/cm²  | 12    | 29.0M       | 6'04''     | 12.5 μs/traj |
|                           | I₀ = 5×10¹⁴ W/cm²  | 12    | 32.9M       | 7'00''     | 12.8 μs/traj |
|                           | I₀ = 7×10¹⁴ W/cm²  | 12    | 35.2M       | 7'35''     | 12.9 μs/traj |
| `test_Molecules`          | H₂ HOMO            | 12    | 14.7M       | 1'45''     | 7.1 μs/traj  |
|                           | CO HOMO            | 12    | 16.2M       | 1'47''     | 6.6 μs/traj  |
|                           | O₂ α-HOMO          | 12    | 15.7M       | 1'45''     | 6.7 μs/traj  |
|                           | O₂ α-HOMO-1        | 12    | 15.3M       | 1'42''     | 6.7 μs/traj  |
|                           | C₆H₆ HOMO          | 12    | 22.3M       | 2'21''     | 6.3 μs/traj  |
|                           | C₆H₆ HOMO-1        | 12    | 22.0M       | 2'19''     | 6.3 μs/traj  |