using HOHQMesh, GLMakie

# Create a new HOHQMesh model project. The project name
# "Mesh_NO1" will be the name of the mesh file
# saved in the directory "NewMeshFiles".
    MY_project = newProject("Mesh_NO1", "Mesh_Files/Meshes_Folder")

# Reset polynomial order of the mesh model curves and output format.
# The "ABAQUS" mesh file format is needed for the adaptive mesh
# capability of Trixi.jl.
    setPolynomialOrder!(MY_project, 3)
    setMeshFileFormat!(MY_project, "ABAQUS")

# add background grid
    addBackgroundGrid!(MY_project, [0.5, 0.5, 0.0])


    line20 = newEndPointsLineCurve("line20", [0,2,0.0], [0,0,0.0])
    line19 = newEndPointsLineCurve("line19", [2,2,0.0], [0,2,0.0])
    line18 = newEndPointsLineCurve("line18", [2,4,0.0], [2,2,0.0])
    line17 = newEndPointsLineCurve("line17", [0,4,0.0], [2,4,0.0])
    line16 = newEndPointsLineCurve("line16", [0,6,0.0], [0,4,0.0])
    line15 = newEndPointsLineCurve("line15", [2,6,0.0], [0,6,0.0])
    line14 = newEndPointsLineCurve("line14", [2,8,0.0], [2,6,0.0])
    line13 = newEndPointsLineCurve("line13", [0,8,0.0], [2,8,0.0])
    line12 = newEndPointsLineCurve("line12", [0,10,0.0], [0,8,0.0])
    line11 = newEndPointsLineCurve("line11", [6,10,0.0], [0,10,0.0])
    line10 = newEndPointsLineCurve("line10", [6,8,0.0], [6,10,0.0])
    line9 = newEndPointsLineCurve("line9", [4,8,0.0], [6,8,0.0])
    line8 = newEndPointsLineCurve("line8", [4,6,0.0], [4,8,0.0])
    line7 = newEndPointsLineCurve("line7", [6,6,0.0], [4,6,0.0])
    line6 = newEndPointsLineCurve("line6", [6,4,0.0], [6,6,0.0])
    line5 = newEndPointsLineCurve("line5", [4,4,0.0], [6,4,0.0])
    line4 = newEndPointsLineCurve("line4", [4,2,0.0], [4,4,0.0])
    line3 = newEndPointsLineCurve("line3", [6,2,0.0], [4,2,0.0])
    line2 = newEndPointsLineCurve("line2", [6,0,0.0], [6,2,0.0])
    line1 = newEndPointsLineCurve("line1", [0,0,0.0], [6,0,0.0])


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
    @info "Model is fully defined."

# add refinement region - postprocessing smoothing - Can comment out to exclude refinement region.
 center = newRefinementCenter("region", "smooth", [3, 5, 0.0], 0.4, 1.0)
 addRefinementRegion!(MY_project, center)

# Visualize the model, refinement region and background grid
# prior to meshing.x

plotProject!(MY_project, MODEL+GRID)

# generate the mesh, and plot the model with grid
# Generate the mesh. Saves the mesh file to the directory "out".
generate_mesh(MY_project)