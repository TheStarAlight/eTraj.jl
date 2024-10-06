using Plots
# using PyPlot
using JLD2
using CodecZlib
using eTraj
using eTraj.Lasers
using Printf
pyplot()

f_ADK   = jldopen("ADK-CTMC_4e14_800nm_cos4_2cyc_CP.jld2")
f_SPANE = jldopen("SPANE-CTMC_4e14_800nm_cos4_2cyc_CP.jld2")
f_SPA   = jldopen("SPA-CTMC_4e14_800nm_cos4_2cyc_CP.jld2")

px = f_ADK["px"]
py = f_ADK["py"]

spec_ADK   = f_ADK["momentum_spec"]
spec_SPANE = f_SPANE["momentum_spec"]
spec_SPA   = f_SPA["momentum_spec"]
prob_ADK   = f_ADK["ion_prob"]
prob_SPANE = f_SPANE["ion_prob"]
prob_SPA   = f_SPA["ion_prob"]
t  = -120:0.5:120
Ax = f_ADK["params"][:laser] |> LaserAx
Ay = f_ADK["params"][:laser] |> LaserAy
Axt = Ax.(t)
Ayt = Ay.(t)
close(f_ADK)
close(f_SPANE)
close(f_SPA)
spec_ADK   .= max.(spec_ADK, 1e-30) / maximum(spec_ADK)
spec_SPANE .= max.(spec_SPANE, 1e-30) / maximum(spec_SPANE)
spec_SPA   .= max.(spec_SPA, 1e-30) / maximum(spec_SPA)

clim = (-3.0,0.0)
xylim = (-2.5,2.5)
fig_x_size = 300
colorbar_x_size = 115
fig_y_size = 300

begin
    fig_ADK = heatmap(px, py, log10.(spec_ADK)', c=:ice)
    plot!(fig_ADK, xlim=xylim, ylim=xylim, clim=clim, size=(fig_x_size,fig_y_size), aspectratio=:equal, tickdirection=:out, colorbar=false)
    plot!(fig_ADK, xlabel=raw"$p_x$ / a.u.", ylabel=raw"$p_y$ / a.u.")
    annotate!(fig_ADK, xylim[1]*0.95,xylim[2]*0.95, text(raw"800↺nm 2-cyc $xy$-CP", :white, :top, :left, 10))
    annotate!(fig_ADK, xylim[1]*0.95,xylim[1]*0.95, text("ADK-CTMC", :white, :bottom, :left, 12))
    annotate!(fig_ADK, xylim[2]*0.95,xylim[1]*0.95, text(@sprintf("rel. prob.\n=%.2f", prob_ADK), :white, :bottom, :right, 9))
    plot!(fig_ADK, -1 .* Axt, -1 .* Ayt, linecolor=RGBA(0.0,0.0,0.0,0.5), legend=false)
    annotate!(fig_ADK, -0.35, 1.2, text(raw"-$\vec{A}(t)$", RGBA(0.0,0.0,0.0,0.5), :bottom, :center, 9))
end

begin
    fig_SPANE = heatmap(px, py, log10.(spec_SPANE)', c=:ice)
    plot!(fig_SPANE, xlim=xylim, ylim=xylim, clim=clim, size=(fig_x_size,fig_y_size), aspectratio=:equal, tickdirection=:out, colorbar=false)
    plot!(fig_SPANE, xlabel=raw"$p_x$ / a.u.", ylabel=raw"$p_y$ / a.u.")
    annotate!(fig_SPANE, xylim[1]*0.95,xylim[2]*0.95, text(raw"800↺nm 2-cyc $xy$-CP", :white, :top, :left, 10))
    annotate!(fig_SPANE, xylim[1]*0.95,xylim[1]*0.95, text("SPANE-CTMC", :white, :bottom, :left, 12))
    annotate!(fig_SPANE, xylim[2]*0.95,xylim[1]*0.95, text(@sprintf("rel. prob.\n=%.2f", prob_SPANE), :white, :bottom, :right, 9))
    plot!(fig_SPANE, -1 .* Axt, -1 .* Ayt, linecolor=RGBA(0.0,0.0,0.0,0.5), legend=false)
    annotate!(fig_SPANE, -0.35, 1.2, text(raw"-$\vec{A}(t)$", RGBA(0.0,0.0,0.0,0.5), :bottom, :center, 9))
end

begin
    fig_SPA = heatmap(px, py, log10.(spec_SPA)', c=:ice)
    plot!(fig_SPA, xlim=xylim, ylim=xylim, clim=clim, size=(fig_x_size,fig_y_size), aspectratio=:equal, tickdirection=:out, colorbar=true)
    plot!(fig_SPA, xlabel=raw"$p_x$ / a.u.", ylabel=raw"$p_y$ / a.u.")
    annotate!(fig_SPA, xylim[1]*0.95,xylim[2]*0.95, text(raw"800↺nm 2-cyc $xy$-CP", :white, :top, :left, 10))
    annotate!(fig_SPA, xylim[1]*0.95,xylim[1]*0.95, text("SPA-CTMC", :white, :bottom, :left, 12))
    annotate!(fig_SPA, xylim[2]*0.95,xylim[1]*0.95, text(@sprintf("rel. prob.\n=%.2f", prob_SPA), :white, :bottom, :right, 9))
    plot!(fig_SPA, -1 .* Axt, -1 .* Ayt, linecolor=RGBA(0.0,0.0,0.0,0.5), legend=false)
    annotate!(fig_SPA, -0.35, 1.2, text(raw"-$\vec{A}(t)$", RGBA(0.0,0.0,0.0,0.5), :bottom, :center, 9))
end

fig = plot(fig_ADK, fig_SPANE, fig_SPA, layout=@layout([a b c]), size=(fig_x_size*3+colorbar_x_size,fig_y_size))
savefig(fig, "figure_2cycs_CP.pdf")
