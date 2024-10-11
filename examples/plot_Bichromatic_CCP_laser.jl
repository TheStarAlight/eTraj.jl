using eTraj
using eTraj.Lasers, eTraj.Units
using Plots
using Colors
# using PyPlot  # PyPlot is implicitly required

l1 = Cos2Laser(peak_int=1e14*W/cm^2, wave_len=800.0nm, cyc_num=8,  ellip= 1.0)
l2 = Cos2Laser(peak_int=1e14*W/cm^2, wave_len=400.0nm, cyc_num=16, ellip=-1.0)
l = BichromaticLaser(l1=l1, l2=l2)

t = -400:0.5:400
Fx = t .|> LaserFx(l)
Fy = t .|> LaserFy(l)
Fmax = max(maximum(Fx), maximum(Fy))
Fx ./= Fmax
Fy ./= Fmax

color = colorant"#009afa"
color_shadow = colorant"#009afa40"
fig = plot(repeat([400],length(t)), Fx, Fy, linecolor=color_shadow)
plot!(fig, t, repeat([1.0],length(Fx)), Fy, linecolor=color_shadow)
plot!(fig, t, Fx, repeat([-1.0],length(Fy)), linecolor=color_shadow)
plot!(fig, t, Fx, Fy, linecolor=color)
plot!(fig, xlim=(-400, 400), ylim=(-1,1), zlim=(-1,1))
plot!(fig, legend=:none, xflip=true, camera=(30,20))
plot!(fig, xlabel=raw"$t$ / a.u.", ylabel=raw"$F_x/F_0$", zlabel=raw"$F_y/F_0$")
plot!(fig, size=(500,400))
savefig(fig, "figure_Bichromatic_laser.pdf")