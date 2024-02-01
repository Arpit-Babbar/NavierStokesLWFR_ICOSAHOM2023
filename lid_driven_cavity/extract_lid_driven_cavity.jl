using Plots, DelimitedFiles
using Trixi, TrixiLW
using TrixiLW: data_dir
elixir = "elixir_navierstokes_lid_driven_cavity_ghia.jl"
import PyPlot as plt
include("../mpl.jl")

sol = trixi_include(joinpath(@__DIR__,elixir),
   time_step_computation = TrixiLW.Adaptive(),
   time_int_tol = 1e-8,
   polydeg = 4)

pd_y = PlotData1D(sol.u[1], semi, slice = :y, point = (0.5,0.5), nvisnodes = 3)

pd_x = PlotData1D(sol.u[1], semi, slice = :x, point = (0.5,0.5), nvisnodes = 3)

# Save extracted data

writedlm("lwfr_x_vs_u2.txt", zip(pd_y.x, pd_y.data[:,2]))

writedlm("lwfr_y_vs_u1.txt", zip(pd_x.x, pd_x.data[:,3]))

# Load and plot extracted data

lw_x = readdlm(joinpath(@__DIR__,"lwfr_x_vs_u2.txt"));

lw_y = readdlm(joinpath(@__DIR__,"lwfr_y_vs_u1.txt"));

pyplot()

gr()

p = plot(lw_x[:,2], lw_x[:,1], width = 3, label = "LW error-based", size = (500,400))
exact_data_y = readdlm(data_dir()*"/lid_driven_cavity_y_vs_u1_r1000.txt")
plot!(p, exact_data_y[:,2], exact_data_y[:,1], seriestype = :scatter, width = 3, label = "Ghia et al.", color = :black)
xlabel!(p, "\$ v_x \$")
ylabel!(p, "\$ y \$")
savefig(p, "y_vs_vx.pdf")
display(p)

p = plot(lw_y[:,2], lw_y[:,1], width = 3, label = "LW error-based", size = (500,400))
exact_data_x = readdlm(data_dir()*"/lid_driven_cavity_x_vs_u2_r1000.txt")
plot!(p, exact_data_x[:,2], exact_data_x[:,1], seriestype = :scatter, color = :black, width = 3,label = "Ghia et al.")
xlabel!(p, "\$ v_y \$")
ylabel!(p, "\$ x \$")
# plot!(flip = true)
savefig(p, "x_vs_vy.pdf")
display(p)

plt.figure()

plt.plot(lw_x[:,2], lw_x[:,1], label = "LWFR")
exact_data_y = readdlm(data_dir()*"/lid_driven_cavity_y_vs_u1_r1000.txt")
plt.scatter(exact_data_y[:,2], exact_data_y[:,1], label = "Ghia et al.", color = "black")
plt.legend()
plt.grid(true)
plt.xlabel("\$ v_x \$")
plt.ylabel("\$ y \$")
plt.savefig("y_vs_vx.pdf")

plt.figure()
plt.plot(lw_y[:,2], lw_y[:,1], label = "LWFR")
plt.grid(true)
exact_data_x = readdlm(data_dir()*"/lid_driven_cavity_x_vs_u2_r1000.txt")
plt.scatter(exact_data_x[:,2], exact_data_x[:,1], label = "Ghia et al.", color = "black")
plt.legend()
plt.xlabel("\$ v_y \$")
plt.ylabel("\$ x \$")
# plot!(flip = true)
plt.savefig("x_vs_vy.pdf")


