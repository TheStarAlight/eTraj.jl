using Plots
# using PyPlot  #! Note: PyPlot is implicitly required, please install PyPlot.jl and configure PyCall.jl properly (see the eTraj doc on how to configure it) before running this script.
using JLD2
using CodecZlib
using Printf
pyplot()

path = [
    "WFAT-CTMC_Hydrogen_HOMO_4e14_800nm_6cyc_CP.jld2",
    "WFAT-CTMC_CarbonMonoxide_HOMO_4e14_800nm_6cyc_CP.jld2",
    "WFAT-CTMC_Oxygen_α-HOMO_4e14_800nm_6cyc_CP.jld2",
    "WFAT-CTMC_Oxygen_α-HOMO-1_4e14_800nm_6cyc_CP.jld2",
    "WFAT-CTMC_Benzene_HOMO_4e14_800nm_6cyc_CP.jld2",
    "WFAT-CTMC_Benzene_HOMO-1_4e14_800nm_6cyc_CP.jld2",
]

orbitals = [
    raw"H₂ HOMO ($1sσ_g$)" * "\n" * raw" (mol. axis $x$)",
    raw"CO HOMO ($3σ_g$)" * "\n" * raw" (mol. axis $x$ (C→O))",
    raw"O₂ α-HOMO ($2pπ_u$)" * "\n" * raw" (nodal planes $xz$ & $yz$)",
    raw"O₂ α-HOMO ($2pπ_u$)" * "\n" * raw" (nodal planes $xy$ & $yz$)",
    raw"C₆H₆ ⌬ HOMO ($2pπ_3$)" * "\n" * raw" (nodal planes $xz$ & $yz$)",
    raw"C₆H₆ ⌬ HOMO ($2pπ_2$)" * "\n" * raw" (nodal planes $xy$ & $yz$)",
]

number = ["(a) ", "(b) ", "(c1) ", "(c2) ", "(d1) ", "(d2) "]

files = jldopen.(path)

px = files[1]["px"]
py = files[1]["py"]
spec = [file["momentum_spec"] for file in files]
prob = [file["ion_prob"] for file in files]
for i in eachindex(spec)
    spec[i] .= max.(spec[i], 1e-30) / maximum(spec[i])
end

clim = (-2.0,0.0)
xylim = (-2.0,2.0)
fig_x_size = 300
colorbar_x_size = 140
fig_y_size = 300

figs = Vector(undef, length(files))

for i in eachindex(files)
    figs[i] = heatmap(px, py, log10.(spec[i])', c=:ice)
    plot!(figs[i], xlim=xylim, ylim=xylim, clim=clim, size=(fig_x_size,fig_y_size), aspectratio=:equal, tickdirection=:out, colorbar=iseven(i))
    plot!(figs[i], xlabel=raw"$p_x$ (a.u.)", ylabel=raw"$p_y$ (a.u.)")
    annotate!(figs[i], xylim[1]*0.95,xylim[2]*0.95, text(number[i] * raw"6-cyc $xy$-CP ↺", :white, :top, :left, 10))
    annotate!(figs[i], xylim[2]*0.95,xylim[2]*0.95, text("WFAT-CTMC", :white, :top, :right, 12))
    annotate!(figs[i], xylim[2]*0.95,xylim[1]*0.95, text(orbitals[i], :white, :bottom, :right, 9))
    annotate!(figs[i], xylim[1]*0.95,xylim[1]*0.95, text(@sprintf("rel. prob.\n=%.1e", prob[i]), :white, :bottom, :left, 9))
end

fig = plot(figs..., layout=@layout([a b; c d; e f]), size=(fig_x_size*2+colorbar_x_size,fig_y_size*3))
savefig(fig, "figure_Molecules.pdf")
