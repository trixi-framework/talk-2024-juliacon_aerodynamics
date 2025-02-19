using OrdinaryDiffEq
using Trixi

###############################################################################
# semidiscretization of the compressible Euler equations

equations = CompressibleEulerEquations2D(1.4)

prandtl_number() = 0.72

Reynolds = 10^7
mu() = 1.4 * 0.8 * 1 / Reynolds

equations_parabolic = CompressibleNavierStokesDiffusion2D(equations, mu = mu(),
                                                          Prandtl = prandtl_number())
                                                     
AoA = 0.02181661564992912 # 1.25 degreee in radians

@inline function initial_condition_mach08_flow(x, t, equations)
    # set the freestream flow parameters
    rho_freestream = 1.4

    #=
    sin_AoA, cos_AoA = sincos(0.02181661564992912)
    v = 0.8

    v1 = cos_AoA * v
    v2 = sin_AoA * v
    =#

    v1 = 0.7998096216639273
    v2 = 0.017451908027648896

    p_freestream = 1.0

    prim = SVector(rho_freestream, v1, v2, p_freestream)
    return prim2cons(prim, equations)
end

initial_condition = initial_condition_mach08_flow

# Boundary conditions for free-stream testing
boundary_condition_free_stream = BoundaryConditionDirichlet(initial_condition)

polydeg = 3

surface_flux = flux_hllc
volume_flux = flux_chandrashekar

basis = LobattoLegendreBasis(polydeg)
shock_indicator = IndicatorHennemannGassner(equations, basis,
                                            alpha_max = 0.5,
                                            alpha_min = 0.001,
                                            alpha_smooth = true,
                                            variable = density_pressure)
volume_integral = VolumeIntegralShockCapturingHG(shock_indicator;
                                                 volume_flux_dg = volume_flux,
                                                 volume_flux_fv = surface_flux)

# DG Solver                                                 
solver = DGSEM(polydeg = polydeg, surface_flux = surface_flux,
               volume_integral = volume_integral)

cd(@__DIR__)
path = "../mesh_data/"
mesh = "NACA4412_2_2D_unique.inp"
mesh_file = path * mesh

boundary_symbols = [:b2_symmetry_y_strong, :b4_farfield_riem, :b5_farfield_riem, 
                    :b7_farfield_riem, :b6_viscous_solid, :b8_to_stitch_a]
mesh = P4estMesh{2}(mesh_file, polydeg = polydeg, boundary_symbols = boundary_symbols)

restart_filename = "../restart_files/restart_viscous_t_1049.h5"

boundary_conditions = Dict(:b2_symmetry_y_strong => boundary_condition_free_stream,
                           :b4_farfield_riem => boundary_condition_free_stream,
                           :b5_farfield_riem => boundary_condition_free_stream,
                           :b7_farfield_riem => boundary_condition_free_stream,
                           :b6_viscous_solid => boundary_condition_slip_wall,
                           :b8_to_stitch_a => boundary_condition_free_stream)

velocity_airfoil = NoSlip((x, t, equations) -> SVector(0.0, 0.0))
heat_airfoil = Adiabatic((x, t, equations) -> 0.0)

boundary_conditions_airfoil = BoundaryConditionNavierStokesWall(velocity_airfoil,
                                                                heat_airfoil)

boundary_conditions_parabolic = Dict(:b2_symmetry_y_strong => boundary_condition_free_stream,
                                     :b4_farfield_riem => boundary_condition_free_stream,
                                     :b5_farfield_riem => boundary_condition_free_stream,
                                     :b7_farfield_riem => boundary_condition_free_stream,
                                     :b6_viscous_solid => boundary_conditions_airfoil,
                                     :b8_to_stitch_a => boundary_condition_free_stream)

semi = SemidiscretizationHyperbolicParabolic(mesh, (equations, equations_parabolic),
                                             initial_condition, solver;
                                             boundary_conditions = (boundary_conditions,
                                                                    boundary_conditions_parabolic))

###############################################################################
# ODE solvers

tspan = (load_time(restart_filename), 10.5)
ode = semidiscretize(semi, tspan, restart_filename)
dt = load_dt(restart_filename)


# Callbacks

summary_callback = SummaryCallback()

analysis_interval = 5000

aoa() = 0.02181661564992912
rho_inf() = 1.4
p_inf() = 1.0
u_inf() = 0.8
l_inf() = 1.0

force_boundary_name = (:b6_viscous_solid,)

include("../src/analysis_pointwise.jl")

pressure_coefficient = AnalysisSurfacePointwise(force_boundary_name,
                                                SurfacePressureCoefficient(p_inf(),
                                                                           rho_inf(),
                                                                           u_inf(),
                                                                           l_inf()))

analysis_callback = Trixi.AnalysisCallback(semi, interval = analysis_interval,
                                           output_directory = "out",
                                           save_analysis = true,
                                           analysis_errors = Symbol[],
                                           analysis_integrals = (),
                                           analysis_pointwise = (pressure_coefficient,))

#=
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     analysis_errors = Symbol[],
                                     analysis_integrals = ())
=#

save_solution = SaveSolutionCallback(interval = analysis_interval,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     solution_variables = cons2prim)

amr_controller = ControllerThreeLevel(semi, shock_indicator,
                                      base_level = 0,
                                      med_level = 1, med_threshold = 0.05,
                                      max_level = 2, max_threshold = 0.1)

amr_callback = AMRCallback(semi, amr_controller,
                           interval = 200,
                           adapt_initial_condition = false)

callbacks = CallbackSet(summary_callback, 
                        save_solution,
                        amr_callback,
                        analysis_callback)

###############################################################################
# run the simulation

sol = solve(ode, SSPRK43(thread = OrdinaryDiffEq.True());
            abstol = 1.0e-8, reltol = 1.0e-8,
            dt = dt,
            ode_default_options()...,
            callback = callbacks);


summary_callback() # print the timer summary
