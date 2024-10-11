using Plots
# using PyPlot  # PyPlot is implicitly required
using JLD2
using CodecZlib
using eTraj
using eTraj.Lasers
pyplot()

f_1e14 = jldopen("ADK-SCTS_Bichromatic_1.0e14_800+400nm_8+16cycs_CounterCP.jld2")
f_3e14 = jldopen("ADK-SCTS_Bichromatic_3.0e14_800+400nm_8+16cycs_CounterCP.jld2")
f_5e14 = jldopen("ADK-SCTS_Bichromatic_5.0e14_800+400nm_8+16cycs_CounterCP.jld2")
f_7e14 = jldopen("ADK-SCTS_Bichromatic_7.0e14_800+400nm_8+16cycs_CounterCP.jld2")

px = f_1e14["px"]
py = f_1e14["py"]
spec_1e14 = f_1e14["momentum_spec"]
spec_3e14 = f_3e14["momentum_spec"]
spec_5e14 = f_5e14["momentum_spec"]
spec_7e14 = f_7e14["momentum_spec"]
t  = -450:1:450
Ax = f_1e14["params"][:laser] |> LaserAx
Ay = f_1e14["params"][:laser] |> LaserAy
A0 = max(maximum(abs, Ax.(t)), maximum(abs, Ay.(t)))
Axt = Ax.(t) / A0
Ayt = Ay.(t) / A0
close(f_1e14)
close(f_3e14)
close(f_5e14)
close(f_7e14)
spec_1e14 .= max.(spec_1e14, 1e-30) / maximum(spec_1e14)
spec_3e14 .= max.(spec_3e14, 1e-30) / maximum(spec_3e14)
spec_5e14 .= max.(spec_5e14, 1e-30) / maximum(spec_5e14)
spec_7e14 .= max.(spec_7e14, 1e-30) / maximum(spec_7e14)

clim = (-3.0,0.0)
xylim = (-2.5,2.5)
fig_x_size = 300
colorbar_x_size = 115
fig_y_size = 300


laser_x_intv = (-2.4,-0.4)
laser_y_intv = (-2.4,-0.4)
scaled_laser_x = @. -Axt*(laser_x_intv[2]-laser_x_intv[1])/2 + (laser_x_intv[2]+laser_x_intv[1])/2  # plot -Ax
scaled_laser_y = @. -Ayt*(laser_y_intv[2]-laser_y_intv[1])/2 + (laser_y_intv[2]+laser_y_intv[1])/2  # plot -Ay

begin
    fig_1e14 = heatmap(px, py, log10.(spec_1e14)', c=:ice)
    plot!(fig_1e14, xlim=xylim, ylim=xylim, clim=clim, size=(fig_x_size,fig_y_size), aspectratio=:equal, tickdirection=:out, colorbar=false)
    plot!(fig_1e14, xlabel=raw"$p_x$ (a.u.)", ylabel=raw"$p_y$ (a.u.)")
    annotate!(fig_1e14, xylim[1]*0.95,xylim[2]*0.95, text(raw"800↺+400↻nm $xy$-CP", :white, :top, :left, 10))
    annotate!(fig_1e14, xylim[2]*0.95,xylim[2]*0.95, text("SCTS", :white, :top, :right, 12))
    annotate!(fig_1e14, xylim[2]*0.95,xylim[1]*0.95, text(raw"$I$₀ = 1.0×10¹⁴ W/cm²", :white, :bottom, :right, 9))
    plot!(fig_1e14, scaled_laser_x, scaled_laser_y, linecolor=:white, legend=false)
    annotate!(fig_1e14, xylim[1]*0.95,xylim[1]*0.5, text(raw"-$\vec{A}(t)$", :white, :bottom, :left, 9))
end
begin
    fig_3e14 = heatmap(px, py, log10.(spec_3e14)', c=:ice)
    plot!(fig_3e14, xlim=xylim, ylim=xylim, clim=clim, size=(fig_x_size+colorbar_x_size,fig_y_size), aspectratio=:equal, tickdirection=:out)
    plot!(fig_3e14, xlabel=raw"$p_x$ (a.u.)", ylabel=raw"$p_y$ (a.u.)")
    annotate!(fig_3e14, xylim[1]*0.95,xylim[2]*0.95, text(raw"800↺+400↻nm $xy$-CP", :white, :top, :left, 10))
    annotate!(fig_3e14, xylim[2]*0.95,xylim[2]*0.95, text("SCTS", :white, :top, :right, 12))
    annotate!(fig_3e14, xylim[2]*0.95,xylim[1]*0.95, text(raw"$I$₀ = 3.0×10¹⁴ W/cm²", :white, :bottom, :right, 9))
end
begin
    fig_5e14 = heatmap(px, py, log10.(spec_5e14)', c=:ice)
    plot!(fig_5e14, xlim=xylim, ylim=xylim, clim=clim, size=(fig_x_size,fig_y_size), aspectratio=:equal, tickdirection=:out, colorbar=false)
    plot!(fig_5e14, xlabel=raw"$p_x$ (a.u.)", ylabel=raw"$p_y$ (a.u.)")
    annotate!(fig_5e14, xylim[1]*0.95,xylim[2]*0.95, text(raw"800↺+400↻nm $xy$-CP", :white, :top, :left, 10))
    annotate!(fig_5e14, xylim[2]*0.95,xylim[2]*0.95, text("SCTS", :white, :top, :right, 12))
    annotate!(fig_5e14, xylim[2]*0.95,xylim[1]*0.95, text(raw"$I$₀ = 5.0×10¹⁴ W/cm²", :white, :bottom, :right, 9))
end
begin
    fig_7e14 = heatmap(px, py, log10.(spec_7e14)', c=:ice)
    plot!(fig_7e14, xlim=xylim, ylim=xylim, clim=clim, size=(fig_x_size+colorbar_x_size,fig_y_size), aspectratio=:equal, tickdirection=:out)
    plot!(fig_7e14, xlabel=raw"$p_x$ (a.u.)", ylabel=raw"$p_y$ (a.u.)")
    annotate!(fig_7e14, xylim[1]*0.95,xylim[2]*0.95, text(raw"800↺+400↻nm $xy$-CP", :white, :top, :left, 10))
    annotate!(fig_7e14, xylim[2]*0.95,xylim[2]*0.95, text("SCTS", :white, :top, :right, 12))
    annotate!(fig_7e14, xylim[2]*0.95,xylim[1]*0.95, text(raw"$I$₀ = 7.0×10¹⁴ W/cm²", :white, :bottom, :right, 9))
end

fig = plot(fig_1e14, fig_3e14, fig_5e14, fig_7e14, layout=@layout([a b; c d]), size=(fig_x_size*2+colorbar_x_size,fig_y_size*2))
savefig(fig, "figure_Bichromatic_CounterCP.pdf")
