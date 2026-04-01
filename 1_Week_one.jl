using OrdinaryDiffEqLowStorageRK
using Trixi
using Pkg
using Plots

###############################################################################
# semidiscretization of the compressible Euler equations

gamma = 1.4
equations = CompressibleEulerEquations2D(gamma)

"""
    initial_condition_sedov_blast_wave(x, t, equations::CompressibleEulerEquations2D)

The Sedov blast wave setup based on example 35.1.4 from Flash
https://flash.rochester.edu/site/flashcode/user_support/flash4_ug_4p8.pdf
"""
function initial_condition_sedov_blast_wave(x, t, equations::CompressibleEulerEquations2D)

    RealT = eltype(x)

    inicenter = SVector(0, 0)
    x_norm = x[1] - inicenter[1]
    y_norm = x[2] - inicenter[2]
    r = sqrt(x_norm^2 + y_norm^2)

    r0 = 0.21875f0
    E = 1

    p0_inner = 3 * (equations.gamma - 1) * E / (3 * convert(RealT, pi) * r0^2)
    p0_outer = convert(RealT, 1.0e-5)

    rho = 1
    v1 = 0
    v2 = 0
    p = r > r0 ? p0_outer : p0_inner

    return prim2cons(SVector(rho, v1, v2, p), equations)
end

initial_condition = initial_condition_sedov_blast_wave

###############################################################################
# Solver setup

surface_flux = FluxLaxFriedrichs(max_abs_speed_naive)
volume_flux = flux_chandrashekar

basis = LobattoLegendreBasis(3)

indicator_sc = IndicatorHennemannGassner(
    equations, basis;
    alpha_max = 0.5,
    alpha_min = 0.001,
    alpha_smooth = true,
    variable = density_pressure
)

volume_integral = VolumeIntegralShockCapturingHG(
    indicator_sc;
    volume_flux_dg = volume_flux,
    volume_flux_fv = surface_flux
)

solver = DGSEM(basis, surface_flux, volume_integral)

coordinates_min = (-2.0, -2.0)
coordinates_max = (2.0, 2.0)

mesh = TreeMesh(
    coordinates_min,
    coordinates_max;
    initial_refinement_level = 6,
    n_cells_max = 100_000,
    periodicity = true
)

semi = SemidiscretizationHyperbolic(
    mesh,
    equations,
    initial_condition,
    solver;
    boundary_conditions = boundary_condition_periodic
)

###############################################################################
# ODE setup

tspan = (0.0, 12.5)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 100
analysis_callback = AnalysisCallback(semi; interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(
    interval = 100,
    save_initial_solution = true,
    save_final_solution = true,
    solution_variables = cons2prim
)

amr_indicator = IndicatorHennemannGassner(
    semi;
    alpha_max = 0.5,
    alpha_min = 0.001,
    alpha_smooth = true,
    variable = density_pressure
)

amr_controller = ControllerThreeLevel(
    semi,
    amr_indicator;
    base_level = 4,
    max_level = 6,
    max_threshold = 0.01
)

amr_callback = AMRCallback(
    semi,
    amr_controller;
    interval = 5,
    adapt_initial_condition = true,
    adapt_initial_condition_only_refine = true
)

stepsize_callback = StepsizeCallback(cfl = 0.8)

###############################################################################
# Pressure probe (currently still using first DOF in memory)

probe_point = (0.0, 0.0)

pressure_history = Float64[]
time_history = Float64[]

function record_pressure!(integrator)

    u = integrator.u

    # WARNING: this is NOT a physical probe
    u_cell = view(u, 1:4)

    prim = cons2prim(u_cell, equations)
    p = prim[4]

    push!(pressure_history, p)
    push!(time_history, integrator.t)
end

condition(u, t, integrator) = true
affect!(integrator) = record_pressure!(integrator)

pressure_callback = DiscreteCallback(condition, affect!)

callbacks = CallbackSet(
    summary_callback,
    analysis_callback,
    alive_callback,
    save_solution,
    amr_callback,
    stepsize_callback,
    pressure_callback
)

###############################################################################
# Run simulation
sol = solve(
    ode,
    CarpenterKennedy2N54(williamson_condition = false);
    dt = 1.0,
    ode_default_options()...,
    callback = callbacks
)

display(plot(sol))

display(plot(
    time_history,
    pressure_history,
    xlabel = "Time",
    ylabel = "Pressure (probe)",
    title = "Probe Pressure vs Time",
    lw = 2
))