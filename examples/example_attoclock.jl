# === example_attoclock.jl ===
#
# With the aim of studying the influence of non-adiabatic effects on the attoclock signal, we employ ADK, SFA-AE and SFA to provide the initial conditions of the electrons, and perform trajectory simulations.
#
# To run the example, a julia intepreter (version ≥ 1.7) with the *SemiclassicalSFI* package installed is required.
# To correctly generate the plots using the code in the example, the HDF5, Colors, ColorSchemes, Plots and PyPlot (which is successfully precompiled) packages are also required to be installed via the package manager of julia.
#
# Usage: place example_attoclock.jl and ijet.jl under the same directory, and execute:
#   julia -t <n_threads> example_attoclock.jl
# where the `n_threads` should be replaced with the desired number of threads allocated for julia.


@info "SemiclassicalSFI Example - Attoclock"

@info "Loading Packages..."
using SemiclassicalSFI
using SemiclassicalSFI.Targets
using SemiclassicalSFI.Lasers

# a Helium atom as target.
t = HeAtom()

l = Cos4Laser(peak_int=4e14, wave_len=800, cyc_num=2, ellip=1.0)

init_cond_list = [:ADK, :SFAAE, :SFA]

for init_cond in init_cond_list
    @info "Running $(init_cond)..."
    filename = "SCSFI_He_4e14_800nm_2cyc_elp1.0_$(init_cond)_CTMC.h5"
    performSFI(
        target = t,
        laser = l,
        init_cond_method = init_cond,
        sample_t_intv = (-80,80),
        sample_t_num = 20000,
        ss_kd_max = 1.5,
        ss_kd_num = 1000,
        ss_kz_max = 1.5,
        ss_kz_num = 50,
        traj_phase_method = :CTMC,
        traj_t_final = 120,
        final_p_max = (2,2,2),
        final_p_num = (500,500,1),
        save_path = filename
    )
end

# ========= PLOTTING =========

@info "Plotting..."
using HDF5
using Plots
pyplot()
include("ijet.jl")
clim = (-5,0)

function plot_heatmap(h5file; kwargs...)
    ion_prob = read(h5file, "ion_prob")
    px = read(h5file, "px")
    py = read(h5file, "py")
    spec = read(h5file, "momentum_spec_2D")
    spec ./= ion_prob
    @. spec = max(1e-100, spec)
    return heatmap(px,py,log10.(spec'), xlim=(-2,2),ylim=(-2,2),aspectratio=:equal, tick_direction=:out,bordercolor=:black,framestyle=:box, c=:ijet,clim=clim ; kwargs...)
end

figure_list = Vector{Any}(nothing,3)
figure_list[1] = plot_heatmap(h5open("SCSFI_He_4e14_800nm_2cyc_elp1.0_ADK_CTMC.h5");    title="CTMC",xlabel=raw"$p_x$ (a.u.)", ylabel=raw"$p_y$ (a.u.)",legend=:none)
figure_list[2] = plot_heatmap(h5open("SCSFI_He_4e14_800nm_2cyc_elp1.0_SFAAE_CTMC.h5");  title="QTMC",xlabel=raw"$p_x$ (a.u.)",legend=:none)
figure_list[3] = plot_heatmap(h5open("SCSFI_He_4e14_800nm_2cyc_elp1.0_SFA_CTMC.h5");    title="SCTS",xlabel=raw"$p_x$ (a.u.)")

l = @layout [a{0.32w} b{0.32w} c{0.36w}]
plot(reshape(figure_list,:)..., size=(900,300), layout=l, fmt=:svg)
savefig("example_attoclock.pdf")
@info "Plot saved as \"example_attoclock.pdf\"."