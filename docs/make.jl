using Documenter, SemiclassicalSFI

makedocs(
        sitename="SemiclassicalSFI.jl",
        pages = [
            "Home"                                              => "index.md",
            "Theory: Initial Conditions"                        => "theory1_initial_conditions.md",
            "Theory: Trajectory Simulation and Phase Methods"   => "theory2_trajectory_simulation_phase_methods.md",
            "Targets"                                           => "program1_targets.md"
        ]
        )