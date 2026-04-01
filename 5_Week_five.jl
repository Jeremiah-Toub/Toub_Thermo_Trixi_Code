###############################################################################
#   initial package calls and setup
    using OrdinaryDiffEqLowStorageRK
    using Trixi
    using Pkg
    using Plots

###############################################################################
#   USER DEFINED  -  Simulation parameters

    # Initial parameters
    Blastwave_center = SVector(3, 1)        # center of the blast wave
    radius_initial = 0.21875                # initial radius of the blast wave
    Energy_initial = 1.0                    # initial energy of the blast wave
    rho = 1                                 # initial density everywhere, including inside the blast wave
    v1 = 0                                  # initial velocity in the x direction everywhere
    v2 = 0                                  # initial velocity in the y direction everywhere
    gamma = 1.4                             # ratio of specific heats for an ideal diatomic gas, which is a common choice for simulating air in fluid dynamics problems. This value is used in the equations of state to relate pressure, density, and internal energy, and it affects the behavior of shock waves and blast waves in the simulation.


    # Resolution and time parameters
    tspan_end = 30                           # end time of the simulation in sim time units
    approximation_degree = 4                # degree of the polynomial approximation used in the DG method. A higher degree allows for a more accurate representation of the solution within each cell, but it also increases the computational cost. Degree of n means that the solution will be approximated using polynomials of degree n within each cell.
    X_resolution = 16                       # number of cells in the x direction for the computational mesh. This determines the spatial resolution of the simulation, with a higher number of cells providing a finer resolution and potentially more accurate results, but also increasing the computational cost.
    Y_resolution = 16                       # number of cells in the y direction for the computational mesh. This determines the spatial resolution of the simulation, with a higher number of cells providing a finer resolution and potentially more accurate results, but also increasing the computational cost.
    interval = 400                          # interval at which to perform analysis and save results during the simulation. This parameter determines how frequently the simulation will output data for analysis and visualization, with a smaller interval providing more frequent updates but also increasing the amount of data generated.     
    numberofgraphs = 16                      # number of graphs to be plotted during the simulation

###############################################################################
#   initializing Function Section

    equations = CompressibleEulerEquations2D(gamma)  # this line initializes the equations of motion for the simulation, specifically the 2D compressible Euler equations with the specified value of gamma. This sets up the mathematical framework for simulating the behavior of the blast wave in the fluid.


    function initial_condition_sedov_blast_wave(x, t, equations::CompressibleEulerEquations2D)  # x is the coordinate vector, t is the time, equations is the PDE system. This is the function that defines the initial condition for the simulation, which will be called at the beginning of the simulation to set up the initial state of the system based on the Sedov blast wave problem.

        RealT = eltype(x)  # get the real type from the input coordinates

        x_norm = x[1] - Blastwave_center[1]  # distance from the center in x direction
        y_norm = x[2] - Blastwave_center[2]  # distance from the center in y direction
        r = sqrt(x_norm^2 + y_norm^2)  # radial distance from the center

        p0_inner = 3 * (equations.gamma - 1) * Energy_initial / (3 * convert(RealT, pi) * radius_initial^2)  # pressure inside the blast wave
        p0_outer = convert(RealT, 1.0e-5)  # pressure outside the blast wave

        #display(typeof(p0_inner))
        #display(typeof(p0_outer))
        p = r > radius_initial ? p0_outer : p0_inner  # if outside the blast wave, the real center point uses p0_outer, otherwise (inside the blast radius) uses p0_inner
        # previous line ensures that the blast wave is initialized with a high pressure inside the radius and a low pressure outside, creating the conditions for the blast wave to propagate outward

        return prim2cons(SVector(rho, v1, v2, p), equations)  # convert the primitive variables (density, velocity, pressure) to conservative variables (density, momentum in x, momentum in y, energy) and return
    end

    initial_condition = initial_condition_sedov_blast_wave  # calls the initial condition function defined above to set up the initial state of the simulation based on the Sedov blast wave problem

###############################################################################
#   Solver setup

    surface_flux = FluxLaxFriedrichs(max_abs_speed_naive)  # this line sets up the numerical flux function for the surface integrals in the DG method. The Lax-Friedrichs flux is a common choice for hyperbolic PDEs, and it uses the maximum absolute wave speed (max_abs_speed_naive) to ensure stability of the numerical scheme. This flux will be used to compute the fluxes across the cell interfaces in the DG method.
    volume_flux = flux_chandrashekar  # this line sets up the numerical flux function for the volume integrals in the DG method

    basis = LobattoLegendreBasis(approximation_degree)  # this line defines the polynomial basis functions used in the DG method. A Lobatto-Legendre basis of degree n correlates to a polynomial degree of n, in this case, degree 3.

    indicator_sc = IndicatorHennemannGassner(  # this line sets up the shock capturing indicator based on the Hennemann-Gassner method. This indicator will be used to identify regions of the solution where shocks or discontinuities are present, and it will help to apply appropriate numerical dissipation in those regions to stabilize the solution. The parameters alpha_max, alpha_min, and alpha_smooth control the behavior of the shock capturing mechanism, while the variable density_pressure indicates that the indicator will be based on a combination of density and pressure to detect shocks.
    equations, basis;
    alpha_max = 0.5,
    alpha_min = 0.001,
    alpha_smooth = false,   # CHECK THIS
    variable = density_pressure
    )


    volume_integral = VolumeIntegralShockCapturingHG(  # this line sets up the volume integral for the DG method, 
    indicator_sc;
    volume_flux_dg = volume_flux,
    volume_flux_fv = surface_flux
    )
    
    # sets up a Discontinuous Galerkin (DG) method that switches to a Finite Volume (FV) method locally where shocks or strong discontinuities are detected by the shock capturing indicator. The volume_flux_dg is used for the DG method in smooth regions, while the volume_flux_fv is used for the FV method in regions where shocks are detected.

    solver = DGSEM(basis, surface_flux, volume_integral)  # this line initializes the solver for the simulation using the Discontinuous Galerkin Spectral Element Method (DGSEM) with the specified basis functions, surface flux, and volume integral. This solver will be used to compute the numerical solution of the compressible Euler equations over time.

###############################################################################
#   Mesh setup

  mesh_file = joinpath("Mesh_Files", "Meshes_Folder/Mesh_NO1.inp")
  mesh = P4estMesh{2}(mesh_file)

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
                                :line20   => boundary_condition_slip_wall
                                )
    #


    semi = SemidiscretizationHyperbolic(
    mesh,
    equations,
    initial_condition,
    solver;
    boundary_conditions = boundary_conditions
    )

###############################################################################
# ODE setup

    tspan = (0.0, tspan_end)
    ode = semidiscretize(semi, tspan)

###############################################################################
# Callbacks

    analysis_interval = interval
    
    summary_callback = SummaryCallback()

    analysis_callback = AnalysisCallback(semi; interval = analysis_interval)
    
    alive_callback = AliveCallback(analysis_interval = analysis_interval)

    save_solution = SaveSolutionCallback(
    interval = interval,
    solution_variables = cons2prim
    )

    stepsize_callback = StepsizeCallback(cfl = 0.8)
    

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

###############################################################################
# Plot results

  # THE FREQUENT ERROR which indicates that "plot is not defined in main" requires you to restart the Julia session and run the code again to properly load the plots package.
  # Easily restart Julia REPL with alt+j then alt+r.

  display(length(sol.u))  # this line displays the total number of time steps taken in the simulation, which can be useful for understanding the computational effort and the resolution of the results.

  datalength = length(sol.u)  # this line stores the total number of time steps taken in the simulation in the variable datalength, which can be used for further analysis or plotting.
  step = PlotData2D(sol.u[10], semi)
  display(plot(step))

  step = PlotData2D(sol.u[100], semi)
  display(plot(step))

  step = PlotData2D(sol.u[200], semi)
  display(plot(step))

  for i in Int(round(datalength/numberofgraphs)):Int(round(datalength/numberofgraphs)):datalength
      step = PlotData2D(sol.u[i], semi)
      display(plot(step))
  end