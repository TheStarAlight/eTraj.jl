using Plots
# using PyPlot  #! Note: PyPlot is implicitly required, please install PyPlot.jl and configure PyCall.jl properly (see the eTraj doc on how to configure it) before running this script.
using JLD2
using CodecZlib
using eTraj
using eTraj.Lasers
pyplot()

fQTMC = jldopen("ADK-QTMC_9e13_800nm_8cyc_LP_ExpRate.jld2")
fSCTS = jldopen("ADK-SCTS_9e13_800nm_8cyc_LP_ExpRate.jld2")

px = fQTMC["px"]
py = fQTMC["py"]
spec_QTMC = fQTMC["momentum_spec"]
spec_SCTS = fSCTS["momentum_spec"]
t  = -450:1:450
Fx = fQTMC["params"][:laser] |> LaserFx
Fxt = Fx.(t)
Fxt ./= maximum(abs, Fxt)
close(fQTMC)
close(fSCTS)
spec_QTMC .= max.(spec_QTMC, 1e-30) / maximum(spec_QTMC)
spec_SCTS .= max.(spec_SCTS, 1e-30) / maximum(spec_SCTS)

clim = (-3.0,0.0)
xylim = (-0.61,0.61)
fig_x_size = 300
colorbar_x_size = 115
fig_y_size = 300

laser_x_intv = (-0.55,-0.1)
laser_y_intv = (-0.55,-0.3)
scaled_laser_x = range(start=laser_x_intv[1], stop=laser_x_intv[2], length=length(t))
scaled_laser_y = @. Fxt*(laser_y_intv[2]-laser_y_intv[1])/2 + (laser_y_intv[2]+laser_y_intv[1])/2

begin
    fig_QTMC = heatmap(px, py, log10.(spec_QTMC), c=:ice)
    plot!(fig_QTMC, xlim=xylim, ylim=xylim, clim=clim, size=(fig_x_size,fig_y_size), aspectratio=:equal, tickdirection=:out, colorbar=false)
    plot!(fig_QTMC, xlabel=raw"$p_y$ (a.u.)", ylabel=raw"$p_x$ (a.u.)")
    annotate!(fig_QTMC, xylim[1]*0.95,xylim[2]*0.95, text(raw"8-cyc $x$-LP", :white, :top, :left, 10))
    annotate!(fig_QTMC, xylim[2]*0.95,xylim[2]*0.95, text("QTMC", :white, :top, :right, 12))
    plot!(fig_QTMC, scaled_laser_x, scaled_laser_y, linecolor=:white, legend=false)
    annotate!(fig_QTMC, xylim[1]*0.95,xylim[1]*0.60, text(raw"$F_x(t)$", :white, :bottom, :left, 9))
end

begin
    fig_SCTS = heatmap(px, py, log10.(spec_SCTS), c=:ice)
    plot!(fig_SCTS, xlim=xylim, ylim=xylim, clim=clim, size=(fig_x_size+colorbar_x_size,fig_y_size), aspectratio=:equal, tickdirection=:out)
    plot!(fig_SCTS, xlabel=raw"$p_y$ (a.u.)", ylabel=raw"$p_x$ (a.u.)")
    annotate!(fig_SCTS, xylim[1]*0.95,xylim[2]*0.95, text(raw"8-cyc $x$-LP", :white, :top, :left, 10))
    annotate!(fig_SCTS, xylim[2]*0.95,xylim[2]*0.95, text("SCTS", :white, :top, :right, 12))
end

fig = plot(fig_QTMC, fig_SCTS, size=(fig_x_size*2+colorbar_x_size,fig_y_size))
savefig(fig, "figure_8cycs_LP.pdf")
