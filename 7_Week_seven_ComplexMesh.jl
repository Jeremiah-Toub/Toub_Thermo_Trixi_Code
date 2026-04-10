###############################################################################
#   initial package calls and setup
    using OrdinaryDiffEqLowStorageRK
    using Trixi
    using Pkg
    using Plots
    using OrdinaryDiffEqSSPRK

###############################################################################
#   USER DEFINED  -  Simulation parameters

    # Initial parameters
    Blastwave_center = SVector(3, 4.5)        # center of the blast wave
    radius_initial = 0.21875                   # initial radius of the blast wave
    Energy_initial = 1.0                       # initial energy of the blast wave
    rho = 1                                    # initial density everywhere
    v1 = 0                                     # initial velocity in the x direction
    v2 = 0                                     # initial velocity in the y direction
    gamma = 1.4                                # ratio of specific heats
    ActiveMesh = "Mesh_NO3_ComplexSquares_TopVent.inp"    # name of the mesh file

#####  WHEN CHANGING MESH ALSO CHANGE THE BOUNDARY CONDITIONS TO MATCH THE CURVE NAMES IN THE MESH FILE. #####

    # Resolution and time parameters
    tspan_end = 3                              # end time of the simulation
    approximation_degree = 3                   # degree of polynomial approximation
    X_resolution = 16                          # number of cells in x direction
    Y_resolution = 16                          # number of cells in y direction
    interval = 400                             # interval for analysis and saving results
    numberofgraphs = 6                         # number of graphs to plot

###############################################################################
#   Initializing Function Section

    equations = CompressibleEulerEquations2D(gamma)

    function initial_condition_sedov_blast_wave(x, t, equations::CompressibleEulerEquations2D)

        RealT = eltype(x)

        x_norm = x[1] - Blastwave_center[1]
        y_norm = x[2] - Blastwave_center[2]
        r = sqrt(x_norm^2 + y_norm^2)

        p0_inner = 3 * (equations.gamma - 1) * Energy_initial / (3 * convert(RealT, pi) * radius_initial^2)
        p0_outer = convert(RealT, 1.0e-5)

        p = r > radius_initial ? p0_outer : p0_inner

        return prim2cons(SVector(rho, v1, v2, p), equations)
    end

    initial_condition = initial_condition_sedov_blast_wave

###############################################################################
#   Solver setup

    surface_flux = FluxLaxFriedrichs(max_abs_speed_naive)
    volume_flux  = flux_chandrashekar

    basis = LobattoLegendreBasis(approximation_degree)

    indicator_sc = IndicatorHennemannGassner(
        equations, basis;
        alpha_max    = 1,
        alpha_min    = 0.001,
        alpha_smooth = true,
        variable     = density_pressure
    )

    volume_integral = VolumeIntegralShockCapturingHG(
        indicator_sc;
        volume_flux_dg = volume_flux,
        volume_flux_fv = surface_flux
    )

solver = DGSEM(basis, surface_flux, volume_integral)
###############################################################################
#   Mesh setup

    mesh_file = joinpath("Mesh_Files", "Meshes_Folder/", ActiveMesh)
    mesh = P4estMesh{2}(mesh_file)

# ============================================================
#  TRIXI BOUNDARY CONDITIONS (reference — paste into sim file)
#
 boundary_conditions = Dict(
     :outer1   => boundary_condition_slip_wall,  # bottom
     :outer2   => boundary_condition_slip_wall,  # right lower
     :outer3   => boundary_condition_slip_wall,  # right notch h
     :outer4   => boundary_condition_slip_wall,  # right notch v
     :outer5   => boundary_condition_slip_wall,  # right notch h
     :outer6   => boundary_condition_slip_wall,  # right upper
     :outer7   => boundary_condition_slip_wall,  # top right
     :vent     => boundary_condition_slip_wall,  # ← replace with inflow/outflow BC
     :outer9   => boundary_condition_slip_wall,  # top left
     :outer10  => boundary_condition_slip_wall,  # left upper
     :outer11  => boundary_condition_slip_wall,  # left notch h
     :outer12  => boundary_condition_slip_wall,  # left notch v
     :outer13  => boundary_condition_slip_wall,  # left notch h
     :outer14  => boundary_condition_slip_wall,  # left lower
     :sq1_b    => boundary_condition_slip_wall,
     :sq1_l    => boundary_condition_slip_wall,
     :sq1_t    => boundary_condition_slip_wall,
     :sq1_r    => boundary_condition_slip_wall,
     :sq2_b    => boundary_condition_slip_wall,
     :sq2_l    => boundary_condition_slip_wall,
     :sq2_t    => boundary_condition_slip_wall,
     :sq2_r    => boundary_condition_slip_wall,
     :sq3_b    => boundary_condition_slip_wall,
     :sq3_l    => boundary_condition_slip_wall,
     :sq3_t    => boundary_condition_slip_wall,
     :sq3_r    => boundary_condition_slip_wall,
     :sq4_b    => boundary_condition_slip_wall,
     :sq4_l    => boundary_condition_slip_wall,
     :sq4_t    => boundary_condition_slip_wall,
     :sq4_r    => boundary_condition_slip_wall,
  )
# ============================================================

    semi = SemidiscretizationHyperbolic(
        mesh,
        equations,
        initial_condition,
        solver;
        boundary_conditions = boundary_conditions
    )

###############################################################################
#   ODE setup

    tspan = (0.0, tspan_end)
    ode = semidiscretize(semi, tspan)

###############################################################################
#   Callbacks

    analysis_interval = interval

    summary_callback  = SummaryCallback()
    analysis_callback = AnalysisCallback(semi; interval = analysis_interval)
    alive_callback    = AliveCallback(analysis_interval = analysis_interval)

    save_solution = SaveSolutionCallback(
        interval           = interval,
        solution_variables = cons2prim
    )

    stepsize_callback = StepsizeCallback(cfl = 0.25)

    callbacks = CallbackSet(
        summary_callback,
        analysis_callback,
        alive_callback,
        save_solution,
        stepsize_callback
    )

###########################f####################################################
#   Run simulation

stage_limiter! = PositivityPreservingLimiterZhangShu(
    thresholds = (5.0e-6, 5.0e-6), 
    variables  = (Trixi.density, pressure)  # note: only density needs Trixi. prefix
)

sol = solve(
    ode,
    CarpenterKennedy2N54(stage_limiter!, williamson_condition = false);  # positional, not keyword
    dt = 1.0e-7,
    ode_default_options()...,
    callback = callbacks,
    saveat = LinRange(0.0, tspan_end, 10)
)

    display(sol.retcode)
    display(sol.t[end])
    display(minimum(sol.u[end]))

###############################################################################
#   Plot results

    # THE FREQUENT ERROR which indicates that "plot is not defined in main" requires
    # you to restart the Julia session and run the code again to properly load the
    # plots package. Easily restart Julia REPL with alt+j then alt+r.

datalength = length(sol.u)
    display(datalength)

# Initial condition snapshot
    plotdata = PlotData2D(sol.u[2], semi)
    display(plot(plotdata))

display("All done!")

indices = round.(Int, LinRange(1, datalength, numberofgraphs))
for i in indices
    local plotdata = PlotData2D(sol.u[i], semi)
    display(plot(plotdata))
end