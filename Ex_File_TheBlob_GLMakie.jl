# Interactive mesh with a curved outer boundary
#
# Create an outer boundary from a set of parametric equations.
# Add manual refinement in a small region around the point (-4, -0.5).
#
# Keywords: outer boundary, parametric equations, refinement center
using Trixi, Plots, GLMakie, Makie

GLMakie.activate!()

using HOHQMesh

# Instantiate the project
blob_project = newProject("TheBlob", "out")

# Create and add the outer boundary curve
xEqn = "x(t) = 4 * cos(2 * pi * t) - 0.6 * cos(8 * pi * t)^3"
yEqn = "y(t) = 4 * sin(2 * pi * t) - 0.5 * sin(11* pi * t)^2"
zEqn = "z(t) = 0.0"
blob = newParametricEquationCurve("Blob", xEqn, yEqn, zEqn)
addCurveToOuterBoundary!(blob_project, blob)

# Add the background grid
addBackgroundGrid!(blob_project, [0.5, 0.5, 0.0])

# Create and add the refinement region
center = newRefinementCenter("region", "smooth", [-4.0, -0.5, 0.0], 0.4, 1.0)
addRefinementRegion!(blob_project, center)

# Generate the mesh

generate_mesh(blob_project)  # creates the mesh but does not display it


plotProject!(blob_project, MODEL+GRID)  # display the model and grid