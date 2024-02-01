using Trixi: trixi_include
using DelimitedFiles

n_levels = 5
min_level = 4
max_level = 3 + n_levels

levels = min_level:max_level

l2_errors = zeros(n_levels)
nxs = zeros(n_levels)

for degree in 3:3
  for (i, level) in enumerate(levels)
    trixi_include(joinpath(@__DIR__,"elixir_advection_diffusion_nonperiodic.jl"),
              initial_refinement_level = level - 1,
      cfl_number = 4.0,
      polydeg = degree)
    l2_errors[i] = analysis_callback(sol).l2[1]
    nxs[i] = 4^(level-1)
  end

  data = hcat(nxs, l2_errors)
  writedlm(joinpath(@__DIR__,"errors_nonperiodic_N$degree.txt"), data)
end

# [3.2323603521345188, 3.48493114277033, 3.5454539982179583]
