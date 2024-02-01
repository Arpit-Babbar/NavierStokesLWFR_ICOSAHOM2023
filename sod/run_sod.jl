using Trixi: trixi_include
using TrixiLW: utils_dir

include(joinpath(utils_dir(), "slice_p4est.jl"))

function run_sod(nx)
  sol = trixi_include("elixir_navierstokes_sod.jl", nx = nx)
  sliced_soln = slice_at_y(sol, 1.0)
  writedlm("soln_nx$nx.txt", sliced_soln)
end

run_sod(70)
run_sod(150)
run_sod(600)

