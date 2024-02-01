using Trixi, TrixiLW
using Plots

###############################################################################
# semidiscretization of the ideal compressible Navier-Stokes equations

# TODO: parabolic; unify names of these accessor functions
prandtl_number() = 0.72
mu() = 2e-5

equations = CompressibleEulerEquations2D(1.4)
equations_parabolic = TrixiLW.CompressibleNavierStokesDiffusion2D(equations, mu=mu(), Prandtl=prandtl_number(),
  gradient_variables=TrixiLW.GradientVariablesConservative())

# Create DG solver with polynomial degree = 3 and (local) Lax-Friedrichs/Rusanov flux as surface flux
solver = DGSEM(polydeg=4, surface_flux=flux_lax_friedrichs,
               volume_integral=TrixiLW.VolumeIntegralFR(TrixiLW.LW()))

coordinates_min = ( 0.0,  -0.5) # minimum coordinates (min(x), min(y))
coordinates_max = ( 1.0,   1.5) # maximum coordinates (max(x), max(y))

# Create a uniformly refined mesh
nx = 70 # 70, 150, 600
trees_per_dimension = (nx, 2)
mesh = P4estMesh(trees_per_dimension,
                 polydeg=4, initial_refinement_level=0,
                 coordinates_min=coordinates_min, coordinates_max=coordinates_max,
                #  periodicity=(false, false),
                 periodicity=(false, true)
                 )

function initial_condition_sod(x, t, equations::CompressibleEulerEquations2D)
  local rho, u, v, p
  if x[1] < 0.5
    rho, u, v, p = 1.0,   0.0, 0.0, 1.0
  else
    rho, u, v, p = 0.125, 0.0, 0.0, 0.1
  end
  return prim2cons(SVector(rho, u, v, p), equations)
end

initial_condition = initial_condition_sod

@inline function boundary_condition_subsonic_constant(U_inner, f_inner, u_inner,
    outer_cache,
    normal_direction::AbstractVector, x, t, dt,
    surface_flux_function, equations::CompressibleEulerEquations2D,
    dg, time_discretization, scaling_factor = 1.0)

    u_boundary = initial_condition_sod(x, t, equations)

    return Trixi.flux_hll(u_inner, u_boundary, normal_direction, equations)
end

# BC types
velocity_bc_lid = NoSlip((x, t, equations) -> SVector(1.0, 0.0))
velocity_bc_cavity = NoSlip((x, t, equations) -> SVector(0.0, 0.0))
heat_bc = Adiabatic((x, t, equations) -> 0.0)
boundary_condition_lid = BoundaryConditionNavierStokesWall(velocity_bc_lid, heat_bc)
boundary_condition_cavity = BoundaryConditionNavierStokesWall(velocity_bc_cavity, heat_bc)

ns_bc_outflow = TrixiLW.OutflowBC((x, t, equations) -> nothing)

# define periodic boundary conditions everywhere
boundary_conditions = (;
                       x_neg = boundary_condition_subsonic_constant,
                       x_pos = boundary_condition_subsonic_constant,
                    #    y_neg = boundary_condition_subsonic_constant,
                    #    y_pos = boundary_condition_subsonic_constant
                      )

boundary_conditions = Dict(
                       :x_neg => boundary_condition_subsonic_constant,
                       :x_pos => boundary_condition_subsonic_constant,
                    #    :y_neg => boundary_condition_subsonic_constant,
                    #    :y_pos => boundary_condition_subsonic_constant
                      )

boundary_conditions_parabolic = Dict( :x_neg => boundary_condition_cavity,
                                    #   :y_neg => boundary_condition_cavity,
                                    #   :y_pos => boundary_condition_lid,
                                      :x_pos => boundary_condition_cavity
                                      )

# A semidiscretization collects data structures and functions for the spatial discretization
semi = TrixiLW.SemidiscretizationHyperbolicParabolic(mesh,
                                                  get_time_discretization(solver),
                                             (equations, equations_parabolic), initial_condition, solver;
                                             boundary_conditions=(boundary_conditions, boundary_conditions_parabolic),
                                             initial_caches = ((;dt = ones(1)), (;)))

###############################################################################
# ODE solvers, callbacks etc.

# Create ODE problem with time span `tspan`
tspan = (0.0, 0.2)
lw_update = TrixiLW.semidiscretize(semi,
                                   get_time_discretization(solver),
                                   tspan);

summary_callback = SummaryCallback()
save_solution = SaveSolutionCallback(interval=1000,
                                     save_initial_solution=true,
                                     save_final_solution=true,
                                     solution_variables=cons2prim)

alive_callback = AliveCallback(alive_interval=100)
analysis_interval = 5000
analysis_callback = AnalysisCallback(semi, interval=analysis_interval)

callbacks = (
  save_solution,
  analysis_callback,
  alive_callback
);

###############################################################################
# run the simulation

time_int_tol = 1e-8
tolerances = (;abstol = time_int_tol, reltol = time_int_tol)
dt_initial = 2.5e-01
cfl_number = 10
sol = TrixiLW.solve_lwfr(lw_update, callbacks, dt_initial, tolerances,
                        #  time_step_computation = TrixiLW.CFLBased(cfl_number),
                         time_step_computation = TrixiLW.Adaptive(),
                        );

return sol;
