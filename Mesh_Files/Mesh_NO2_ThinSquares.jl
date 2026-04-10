using HOHQMesh, GLMakie

# Create a new HOHQMesh model project. The project name
# "Mesh_NO1" will be the name of the mesh file
# saved in the directory "NewMeshFiles".
    MY_project = newProject("Mesh_NO2_ThinSquares", "Mesh_Files/Meshes_Folder")

# Reset polynomial order of the mesh model curves and output format.
# The "ABAQUS" mesh file format is needed for the adaptive mesh
# capability of Trixi.jl.
    setPolynomialOrder!(MY_project, 3)
    setMeshFileFormat!(MY_project, "ABAQUS")

# add background grid
    addBackgroundGrid!(MY_project, [0.25, 0.25, 0.0])


line1 = newEndPointsLineCurve("line1", [0,0,0.0], [3,0,0.0])
line2 = newEndPointsLineCurve("line2", [3,0,0.0], [3,1.25,0.0])
line3 = newEndPointsLineCurve("line3", [3,1.25,0.0], [2.5,1.25,0.0])
line4 = newEndPointsLineCurve("line4", [2.5,1.25,0.0], [2.5,1.75,0.0])
line5 = newEndPointsLineCurve("line5", [2.5,1.75,0.0], [3,1.75,0.0])
line6 = newEndPointsLineCurve("line6", [3,1.75,0.0], [3,3.25,0.0])
line7 = newEndPointsLineCurve("line7", [3,3.25,0.0], [2.5,3.25,0.0])
line8 = newEndPointsLineCurve("line8", [2.5,3.25,0.0], [2.5,3.75,0.0])
line9 = newEndPointsLineCurve("line9", [2.5,3.75,0.0], [3,3.75,0.0])
line10 = newEndPointsLineCurve("line10", [3,3.75,0.0], [3,5,0.0])
line11 = newEndPointsLineCurve("line11", [3,5,0.0], [0,5,0.0])
line12 = newEndPointsLineCurve("line12", [0,5,0.0], [0,3.75,0.0])
line13 = newEndPointsLineCurve("line13", [0,3.75,0.0], [0.5,3.75,0.0])
line14 = newEndPointsLineCurve("line14", [0.5,3.75,0.0], [0.5,3.25,0.0])
line15 = newEndPointsLineCurve("line15", [0.5,3.25,0.0], [0,3.25,0.0])
line16 = newEndPointsLineCurve("line16", [0,3.25,0.0], [0,1.75,0.0])
line17 = newEndPointsLineCurve("line17", [0,1.75,0.0], [0.5,1.75,0.0])
line18 = newEndPointsLineCurve("line18", [0.5,1.75,0.0], [0.5,1.25,0.0])
line19 = newEndPointsLineCurve("line19", [0.5,1.25,0.0], [0,1.25,0.0])
line20 = newEndPointsLineCurve("line20", [0,1.25,0.0], [0,0,0.0])


line21 = newEndPointsLineCurve("line21", [1.0,2.0,0.0], [2.0,2.0,0.0])   # bottom edge
line22 = newEndPointsLineCurve("line22", [2.0,2.0,0.0], [2.0,3.0,0.0])   # right edge
line23 = newEndPointsLineCurve("line23", [2.0,3.0,0.0], [1.0,3.0,0.0])   # top edge
line24 = newEndPointsLineCurve("line24", [1.0,3.0,0.0], [1.0,2.0,0.0])   # left edge



    addCurveToOuterBoundary!(MY_project, line1)
    addCurveToOuterBoundary!(MY_project, line2)
    addCurveToOuterBoundary!(MY_project, line3)
    addCurveToOuterBoundary!(MY_project, line4)
    addCurveToOuterBoundary!(MY_project, line5)
    addCurveToOuterBoundary!(MY_project, line6)
    addCurveToOuterBoundary!(MY_project, line7)
    addCurveToOuterBoundary!(MY_project, line8)
    addCurveToOuterBoundary!(MY_project, line9)
    addCurveToOuterBoundary!(MY_project, line10)
    addCurveToOuterBoundary!(MY_project, line11)
    addCurveToOuterBoundary!(MY_project, line12)
    addCurveToOuterBoundary!(MY_project, line13)
    addCurveToOuterBoundary!(MY_project, line14)
    addCurveToOuterBoundary!(MY_project, line15)
    addCurveToOuterBoundary!(MY_project, line16)
    addCurveToOuterBoundary!(MY_project, line17)
    addCurveToOuterBoundary!(MY_project, line18)
    addCurveToOuterBoundary!(MY_project, line19)
    addCurveToOuterBoundary!(MY_project, line20)

    addCurveToInnerBoundary!(MY_project, line21, "squareLoop")
    addCurveToInnerBoundary!(MY_project, line22, "squareLoop")
    addCurveToInnerBoundary!(MY_project, line23, "squareLoop")
    addCurveToInnerBoundary!(MY_project, line24, "squareLoop")
    @info "Model is fully defined."

# add refinement region - postprocessing smoothing - Can comment out to exclude refinement region.
 center = newRefinementCenter("region", "smooth", [1.5, 2.5, 0.0], 0.4, 1.0)
 addRefinementRegion!(MY_project, center)

# Visualize the model, refinement region and background grid
# prior to meshing.x

plotProject!(MY_project, MODEL+GRID)

# generate the mesh, and plot the model with grid
# Generate the mesh. Saves the mesh file to the directory "out".
generate_mesh(MY_project)


#=
boundary conditions to match the curve names in the mesh file. #####
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
                                =#