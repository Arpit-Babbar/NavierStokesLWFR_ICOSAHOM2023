using Trixi: trixi_include
using DelimitedFiles

n_levels = 4
min_level = 4
max_level = 3 + n_levels

levels = min_level:max_level

l2_errors = zeros(n_levels)
nxs = zeros(n_levels)

for degree in (2,4)
  for (i, level) in enumerate(levels)
    trixi_include(joinpath(@__DIR__,"elixir_navierstokes_convergence.jl"),
              initial_refinement_level = level-1,
      cfl_number = 4.0,
      polydeg = degree)
    l2_errors[i] = analysis_callback(sol).l2[1]
    nxs[i] = 4^(level-1)
  end

  data = hcat(nxs, l2_errors)
  writedlm(joinpath(@__DIR__,"errors_NS_N$degree.txt"), data)
end
