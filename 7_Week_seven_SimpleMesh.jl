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
    Blastwave_center = SVector(0.5, .75)        # center of the blast wave
    radius_initial = 0.21875                   # initial radius of the blast wave
    Energy_initial = 3.0                       # initial energy of the blast wave
    rho = 1                                    # initial density everywhere
    v1 = 0                                     # initial velocity in the x direction
    v2 = 0                                     # initial velocity in the y direction
    gamma = 1.4                                # ratio of specific heats
    ActiveMesh = "Mesh_NO2_ThinSquares.inp"    # name of the mesh file

#####  WHEN CHANGING MESH ALSO CHANGE THE BOUNDARY CONDITIONS TO MATCH THE CURVE NAMES IN THE MESH FILE. #####

    # Resolution and time parameters
    tspan_end = 13                              # end time of the simulation
    approximation_degree = 3                   # degree of polynomial approximation
    X_resolution = 16                          # number of cells in x direction
    Y_resolution = 16                          # number of cells in y direction
    interval = 400                             # interval for analysis and saving results
    numberofgraphs = 15                         # number of graphs to plot
#   early_snapshot = 1e-6                       # a few timesteps in (adjust as needed)

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
 boundary_conditions = Dict( :line1    => boundary_condition_slip_wall,
        :line2    => boundary_condition_slip_wall,
        :line3    => boundary_condition_slip_wall,
        :line4    => boundary_condition_slip_wall,
        :line5    => boundary_condition_slip_wall,
        :line6    => boundary_condition_slip_wall,
        :line7    => boundary_condition_slip_wall,
        :line8    => boundary_condition_slip_wall,
        :line9    => boundary_condition_slip_wall,
        :line10   => boundary_condition_slip_wall,
        :line11   => boundary_condition_slip_wall,
        :line12   => boundary_condition_slip_wall,
        :line13   => boundary_condition_slip_wall,
        :line14   => boundary_condition_slip_wall,
        :line15   => boundary_condition_slip_wall,                           
        :line16   => boundary_condition_slip_wall,
        :line17   => boundary_condition_slip_wall,                    
        :line18   => boundary_condition_slip_wall,
        :line19   => boundary_condition_slip_wall,
        :line20   => boundary_condition_slip_wall,
        :line21   => boundary_condition_slip_wall,
        :line22   => boundary_condition_slip_wall,
        :line23   => boundary_condition_slip_wall,
        :line24   => boundary_condition_slip_wall
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
    # Previously at cfl = 0.5

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

    #Spacedgraphs = numberofgraphs-1
    #saveat_times = sort(vcat(early_snapshot, LinRange(0.0, tspan_end, Spacedgraphs)))
    saveat_times = LinRange(0.0, tspan_end, numberofgraphs)

    sol = solve(
        ode,
        CarpenterKennedy2N54(stage_limiter!, williamson_condition = false);
        dt = 1.0e-7,
        ode_default_options()...,
        callback = callbacks,
        saveat = saveat_times
    )

###############################################################################
#   Troubleshooting outputs
    display(sol.retcode)
    display(sol.t[end])
    display(minimum(sol.u[end]))

    datalength = length(sol.u)
    display(datalength)

###############################################################################
#   Plot results


    # THE FREQUENT ERROR which indicates that "plot is not defined in main" requires
    # you to restart the Julia session and run the code again to properly load the
    # plots package. Easily restart Julia REPL with alt+j then alt+r.


    function get_p_values(u, semi)
        pds = PlotData2D(u, semi)["p"]
        vid = pds.variable_id
        return StructArrays.component(pds.plot_data.data, vid)
    end

    p_min = minimum(minimum(get_p_values(u, semi)) for u in sol.u[2:end])
    p_max = maximum(maximum(get_p_values(u, semi)) for u in sol.u[2:end])

    # Early snapshot with its own scale
    display(plot(PlotData2D(sol.u[1], semi)["p"], title = "Pressure t=$(round(sol.t[1], sigdigits=3))", color = :turbo))

    # Rest with shared scale
    for i in 2:length(sol.u)
        local plotdata = PlotData2D(sol.u[i], semi)
        display(plot(plotdata["p"], clims = (p_min, p_max), 
            title = "Pressure t=$(round(sol.t[i], sigdigits=3))",
            color = :turbo))
    end