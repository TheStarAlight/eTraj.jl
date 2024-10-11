using eTraj
using eTraj.Targets
using Plots
# using PyPlot  # PyPlot is implicitly required
pyplot()

molH2 = get_mol("Hydrogen")
molCO = get_mol("Carbon Monoxide")
molO2 = get_mol("Oxygen")
molC6H6 = get_mol("Benzene")

θ_range = 0.0:π/90:π
χ_range = repeat([π/2], length(θ_range))
fig_x_size = 400
fig_y_size = 300

begin
    H2_G00sq = MolWFATStructureFactor_G(molH2, 0, 0, 0, θ_range, χ_range) .|> abs2
    figH2 = plot(θ_range .* (180/π), H2_G00sq)
    plot!(figH2, xlabel=raw"$θ$ / deg", ylabel=raw"HOMO $|G_{00}|^2$ / a.u.")
    plot!(figH2, legend=:none, size=(fig_x_size, fig_y_size), framestyle=:box)
    plot!(figH2, xticks=0:45:180, xlim=(0,180), ylim=(0,4))
    annotate!(figH2, 0.97*180, 0.03*4, text("H₂", :black, :bottom, :right, 24))
    annotate!(figH2, 0.03*180, 0.97*4, text("(a)", :black, :top, :left, 12))
end

begin
    CO_G00sq = MolWFATStructureFactor_G(molCO, 0, 0, 0, θ_range, χ_range) .|> abs2
    figCO = plot(θ_range .* (180/π), CO_G00sq)
    plot!(figCO, xlabel=raw"$θ$ / deg", ylabel=raw"HOMO $|G_{00}|^2$ / a.u.")
    plot!(figCO, legend=:none, size=(fig_x_size, fig_y_size), framestyle=:box)
    plot!(figCO, xticks=0:45:180, xlim=(0,180), ylim=(0,70))
    annotate!(figCO, 0.97*180, 0.03*70, text("CO", :black, :bottom, :right, 24))
    annotate!(figCO, 0.03*180, 0.97*70, text("(b)", :black, :top, :left, 12))
end

begin
    O2HOMO_G00   = MolWFATStructureFactor_G(molO2, (1,0) , 0, 0, θ_range, χ_range) .|> abs2
    O2HOMOm1_G01 = MolWFATStructureFactor_G(molO2, (1,-1), 0, 1, θ_range, χ_range) .|> abs2
    figO2 = plot(θ_range .* (180/π), O2HOMO_G00,   label=raw"HOMO   , $ν$=(0,0)")
    plot!(figO2, θ_range .* (180/π), O2HOMOm1_G01, label=raw"HOMO-1, $ν$=(0,1)")
    plot!(figO2, xlabel=raw"$θ$ / deg", ylabel=raw"$|G_ν|^2$ / a.u.")
    plot!(figO2, legend=:top, size=(fig_x_size, fig_y_size), framestyle=:box)
    plot!(figO2, xticks=0:45:180, xlim=(0,180), ylim=(0,12))
    annotate!(figO2, 0.97*180, 0.97*12, text("O₂", :black, :top, :right, 24))
    annotate!(figO2, 0.03*180, 0.97*12, text("(c) " * raw"$χ$=90°", :black, :top, :left, 12))
end

begin
    C6H6HOMO_G00 = MolWFATStructureFactor_G(molC6H6, 0, 0, 0, θ_range, χ_range) .|> abs2
    C6H6HOMOm1_G01 = MolWFATStructureFactor_G(molC6H6, -1, 0, 1, θ_range, χ_range) .|> abs2
    figC6H6 = plot(θ_range .* (180/π), C6H6HOMO_G00,  label=raw"HOMO   , $ν$=(0,0)")
    plot!(figC6H6, θ_range .* (180/π), C6H6HOMOm1_G01, label=raw"HOMO-1, $ν$=(0,1)")
    plot!(figC6H6, xlabel=raw"$θ$ / deg", ylabel=raw"$|G_ν|^2$ / a.u.")
    plot!(figC6H6, legend=:top, size=(fig_x_size, fig_y_size), framestyle=:box)
    plot!(figC6H6, xticks=0:45:180, xlim=(0,180), ylim=(0,25))
    annotate!(figC6H6, 0.97*180, 0.97*25, text("C₆H₆", :black, :top, :right, 24))
    annotate!(figC6H6, 0.03*180, 0.97*25, text("(d) " * raw"$χ$=90°", :black, :top, :left, 12))
end

fig = plot(figH2, figCO, figO2, figC6H6, layout=@layout([a b;c d]), size=(2*fig_x_size, 2*fig_y_size))
savefig(fig, "figure_Molecules_G.pdf")
