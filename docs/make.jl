using Documenter, SemiclassicalSFI

makedocs(
        sitename="SemiclassicalSFI.jl",
        pages = [
            "Home"      => "index.md",
            "Theory"    => [
                "Initial Conditions"                        => "theory1_initial_conditions.md",
                "Trajectory Simulation and Phase Methods"   => "theory2_trajectory_simulation_phase_methods.md"
                ],
            "Manual"    => [
                "Targets"   => "manual1_targets.md"
                ]
        ],
        doctest = false
        )