using Documenter, eTraj, eTraj.Lasers, eTraj.Targets, eTraj.Units

makedocs(
    sitename = "eTraj.jl",
    repo = Remotes.GitHub("TheStarAlight", "eTraj.jl"),
    pages = [
        "Home"      => "index.md",
        "Step-by-step Tutorial" => "stepbystep_tutorial.md",
        "Troubleshooting" => "troubleshooting.md",
        "Theory"    => [
            "Initial Conditions"                        => "theory1_init_cond.md",
            "Trajectory Simulation and Phase Methods"   => "theory2_traj_simu_quant_phase.md"
            ],
        "Manual"    => [
            "The `Lasers` Module"        => "manual1_lasers.md",
            "The `Targets` Module"       => "manual2_targets.md",
            "Sampling & Traj Simulation" => "manual3_sampling_traj_simu.md"
            ],
        "Examples"  => [
            "Overview"                              => "example_overview.md",
            "Attoclock & Initial Condition Methods" => "example1_attoclock.md",
            "Short LP Pulses & Phase Methods"       => "example2_ShortLP_phase_methods.md",
            "Bichromatic ω-2ω Pulses"               => "example3_bichromatic_pulse.md",
            "WFAT-CTMC for Molecules"               => "example4_WFAT-CTMC.md"
        ]
    ],
    format = Documenter.HTML(
        ansicolor=true,
        edit_link=nothing,
        footer="· *eTraj.jl* Documentation · by *Mingyu Zhu* and other contributors",
        mathengine=
        Documenter.KaTeX(
                Dict(:macros =>
                    Dict(
                        raw"\rm"    => raw"\mathrm",
                        raw"\pd"    => raw"\partial",
                        raw"\dd"    => raw"\mathrm{d}",
                        raw"\ii"    => raw"\mathrm{i}",
                        raw"\ee"    => raw"\mathrm{e}",
                        raw"\Re"    => raw"\mathcal{R}",
                        raw"\Im"    => raw"\mathcal{I}",
                        raw"\rr"    => raw"\bm{r}",
                        raw"\pp"    => raw"\bm{p}",
                        raw"\kk"    => raw"\bm{k}",
                        raw"\aa"    => raw"\bm{a}",
                        raw"\AA"    => raw"\bm{A}",
                        raw"\FF"    => raw"\bm{F}",
                        raw"\LL"    => raw"\bm{L}",
                        raw"\RRh"   => raw"\hat{\bm{R}}",
                        raw"\Ip"    => raw"I_{\mathrm{p}}",
                        raw"\Sp"    => raw"S_{\pp}",
                        raw"\ti"    => raw"t_{\mathrm{i}}",
                        raw"\tf"    => raw"t_{\mathrm{f}}",
                        raw"\ts"    => raw"t_{\mathrm{s}}",
                        raw"\tr"    => raw"t_{\mathrm{r}}",
                        raw"\kp"    => raw"k_{\perp}",
                        raw"\kt"    => raw"k_{\mathrm{t}}",
                        raw"\kkt"   => raw"\kk_{\mathrm{t}}",
                        # KaTeX special
                        raw"\mel"   => raw"\braket{#1 | #2 | #3}", # KaTeX doesn't support \mel, use \braket instead
                        raw"\abs"   => raw"\lvert #1 \rvert",
                        raw"\bk"    => raw"\braket{#1 | #2}",
                    ),
                )
            )
        )
    )