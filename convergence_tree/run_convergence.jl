using Trixi: trixi_include
using DelimitedFiles

n_levels = 5
min_level = 4
max_level = 3 + n_levels

levels = min_level:max_level

trees_per_dimension = (8,8)

l2_errors = zeros(n_levels)
nxs = zeros(n_levels)

println("Running elixir_advection_diffusion_1em6")
for degree in 2:2
  for (i, level) in enumerate(levels)
    trixi_include(joinpath(@__DIR__,"elixir_advection_diffusion_1em6.jl"),
              initial_refinement_level = level-1,
      cfl_number = 0.98,
      polydeg = degree)
    l2_errors[i] = analysis_callback(sol).l2[1]
    nxs[i] = 4^(level-1)
  end

  data = hcat(nxs, l2_errors)
  writedlm(joinpath(@__DIR__,"errors_1em6_N$degree.txt"), data)
end

n_levels = 5
min_level = 4
max_level = 3 + n_levels

levels = min_level:max_level

l2_errors = zeros(n_levels)
nxs = zeros(n_levels)
println("Running elixir_advection_diffusion_5em2")
for degree in 2:2
  for (i, level) in enumerate(levels)
    trixi_include(joinpath(@__DIR__,"elixir_advection_diffusion_5em2.jl"),
              initial_refinement_level = level-1,
      cfl_number = 5.0, # cfl_number = 4 gives 4.94 convergence rate
      tspan = (0.0, 1.0),
      polydeg = degree)
      l2_errors[i] = analysis_callback(sol).l2[1]
      nxs[i] = 4^(level-1)
  end
  data = hcat(nxs, l2_errors)
  writedlm(joinpath(@__DIR__,"errors_5em2_N$degree.txt"), data)
end

println("Running elixir_advection_diffusion_1em12")
for degree in 2:2
  for (i, level) in enumerate(levels)
    trixi_include(joinpath(@__DIR__,"elixir_advection_diffusion_1em12.jl"),
                  initial_refinement_level = level-1,
                  cfl_number = 0.98,
                  polydeg = degree)
      l2_errors[i] = analysis_callback(sol).l2[1]
      nxs[i] = 4^(level-1)
  end
  data = hcat(nxs, l2_errors)
  writedlm(joinpath(@__DIR__,"errors_1em12_N$degree.txt"), data)
end
