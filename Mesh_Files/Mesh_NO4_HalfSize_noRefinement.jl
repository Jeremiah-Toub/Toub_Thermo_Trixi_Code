using HOHQMesh, GLMakie

# ============================================================
#  Project setup
# ============================================================
MY_project = newProject("Mesh_NO4_HalfSize_noRefinement", "Mesh_Files/Meshes_Folder")

setPolynomialOrder!(MY_project, 2)
setMeshFileFormat!(MY_project, "ABAQUS")

addBackgroundGrid!(MY_project, [.125, .125, 0.0])   # was [0.25, 0.25, 0.0]

# ============================================================
#  OUTER BOUNDARY
#
#  Domain is a 3 × 4.5 rectangle (x: 0→3, y: 0→4.5) with:
#    - A notch cut into the LEFT  wall at y = 1.5 → 2.25
#    - A notch cut into the RIGHT wall at y = 1.5 → 2.25
#    - A vent OPENING at the TOP  from x = 1.25 → 1.75
#      (open boundary — no slip wall on that segment)
#
#  All coordinates are exactly half of the original.
#  Traversal order: counter-clockwise starting at origin.
# ============================================================

# --- Bottom edge ---
outer1  = newEndPointsLineCurve("outer1",  [0.0, 0.0, 0.0], [3.0, 0.0, 0.0])

# --- Right wall: bottom section ---
outer2  = newEndPointsLineCurve("outer2",  [3.0, 0.0, 0.0], [3.0, 1.5, 0.0])

# --- Right notch (indents inward, i.e. toward x = 2.75) ---
outer3  = newEndPointsLineCurve("outer3",  [3.0, 1.5,  0.0], [2.75, 1.5,  0.0])
outer4  = newEndPointsLineCurve("outer4",  [2.75, 1.5, 0.0], [2.75, 2.25, 0.0])
outer5  = newEndPointsLineCurve("outer5",  [2.75, 2.25, 0.0], [3.0, 2.25, 0.0])

# --- Right wall: top section ---
outer6  = newEndPointsLineCurve("outer6",  [3.0, 2.25, 0.0], [3.0, 4.5, 0.0])

# --- Top edge: right portion (x = 3 → 1.75) ---
outer7  = newEndPointsLineCurve("outer7",  [3.0, 4.5, 0.0], [1.75, 4.5, 0.0])

# --- Vent opening at top (x = 1.75 → 1.25): open boundary ---
outer8  = newEndPointsLineCurve("vent",    [1.75, 4.5, 0.0], [1.25, 4.5, 0.0])

# --- Top edge: left portion (x = 1.25 → 0) ---
outer9  = newEndPointsLineCurve("outer9",  [1.25, 4.5, 0.0], [0.0, 4.5, 0.0])

# --- Left wall: top section ---
outer10 = newEndPointsLineCurve("outer10", [0.0, 4.5, 0.0], [0.0, 2.25, 0.0])

# --- Left notch (indents inward, i.e. toward x = 0.25) ---
outer11 = newEndPointsLineCurve("outer11", [0.0,  2.25, 0.0], [0.25, 2.25, 0.0])
outer12 = newEndPointsLineCurve("outer12", [0.25, 2.25, 0.0], [0.25, 1.5,  0.0])
outer13 = newEndPointsLineCurve("outer13", [0.25, 1.5,  0.0], [0.0,  1.5,  0.0])

# --- Left wall: bottom section ---
outer14 = newEndPointsLineCurve("outer14", [0.0, 1.5, 0.0], [0.0, 0.0, 0.0])

# Add all outer curves
addCurveToOuterBoundary!(MY_project, outer1)
addCurveToOuterBoundary!(MY_project, outer2)
addCurveToOuterBoundary!(MY_project, outer3)
addCurveToOuterBoundary!(MY_project, outer4)
addCurveToOuterBoundary!(MY_project, outer5)
addCurveToOuterBoundary!(MY_project, outer6)
addCurveToOuterBoundary!(MY_project, outer7)
addCurveToOuterBoundary!(MY_project, outer8)
addCurveToOuterBoundary!(MY_project, outer9)
addCurveToOuterBoundary!(MY_project, outer10)
addCurveToOuterBoundary!(MY_project, outer11)
addCurveToOuterBoundary!(MY_project, outer12)
addCurveToOuterBoundary!(MY_project, outer13)
addCurveToOuterBoundary!(MY_project, outer14)

@info "Outer boundary defined."

# ============================================================
#  INNER BOUNDARIES — four square holes
#
#  Each hole is traversed CLOCKWISE (required by HOHQMesh for
#  inner boundaries so the fluid domain stays on the left).
#
#  Original centres and half-widths halved:
#    sq1: centre (0.75, 3.5),  half = 0.3  → x: 0.45→1.05, y: 3.2→3.8
#    sq2: centre (2.25, 3.5),  half = 0.3  → x: 1.95→2.55, y: 3.2→3.8
#    sq3: centre (0.75, 1.0),  half = 0.3  → x: 0.45→1.05, y: 0.7→1.3
#    sq4: centre (2.25, 1.0),  half = 0.3  → x: 1.95→2.55, y: 0.7→1.3
# ============================================================

# --- Square hole 1 (top-left) ---
sq1_b = newEndPointsLineCurve("sq1_b", [1.05, 3.2, 0.0], [0.45, 3.2, 0.0])
sq1_l = newEndPointsLineCurve("sq1_l", [0.45, 3.2, 0.0], [0.45, 3.8, 0.0])
sq1_t = newEndPointsLineCurve("sq1_t", [0.45, 3.8, 0.0], [1.05, 3.8, 0.0])
sq1_r = newEndPointsLineCurve("sq1_r", [1.05, 3.8, 0.0], [1.05, 3.2, 0.0])

addCurveToInnerBoundary!(MY_project, sq1_b, "squareLoop1")
addCurveToInnerBoundary!(MY_project, sq1_l, "squareLoop1")
addCurveToInnerBoundary!(MY_project, sq1_t, "squareLoop1")
addCurveToInnerBoundary!(MY_project, sq1_r, "squareLoop1")

# --- Square hole 2 (top-right) ---
sq2_b = newEndPointsLineCurve("sq2_b", [2.55, 3.2, 0.0], [1.95, 3.2, 0.0])
sq2_l = newEndPointsLineCurve("sq2_l", [1.95, 3.2, 0.0], [1.95, 3.8, 0.0])
sq2_t = newEndPointsLineCurve("sq2_t", [1.95, 3.8, 0.0], [2.55, 3.8, 0.0])
sq2_r = newEndPointsLineCurve("sq2_r", [2.55, 3.8, 0.0], [2.55, 3.2, 0.0])

addCurveToInnerBoundary!(MY_project, sq2_b, "squareLoop2")
addCurveToInnerBoundary!(MY_project, sq2_l, "squareLoop2")
addCurveToInnerBoundary!(MY_project, sq2_t, "squareLoop2")
addCurveToInnerBoundary!(MY_project, sq2_r, "squareLoop2")

# --- Square hole 3 (bottom-left) ---
sq3_b = newEndPointsLineCurve("sq3_b", [1.05, 0.7, 0.0], [0.45, 0.7, 0.0])
sq3_l = newEndPointsLineCurve("sq3_l", [0.45, 0.7, 0.0], [0.45, 1.3, 0.0])
sq3_t = newEndPointsLineCurve("sq3_t", [0.45, 1.3, 0.0], [1.05, 1.3, 0.0])
sq3_r = newEndPointsLineCurve("sq3_r", [1.05, 1.3, 0.0], [1.05, 0.7, 0.0])

addCurveToInnerBoundary!(MY_project, sq3_b, "squareLoop3")
addCurveToInnerBoundary!(MY_project, sq3_l, "squareLoop3")
addCurveToInnerBoundary!(MY_project, sq3_t, "squareLoop3")
addCurveToInnerBoundary!(MY_project, sq3_r, "squareLoop3")

# --- Square hole 4 (bottom-right) ---
sq4_b = newEndPointsLineCurve("sq4_b", [2.55, 0.7, 0.0], [1.95, 0.7, 0.0])
sq4_l = newEndPointsLineCurve("sq4_l", [1.95, 0.7, 0.0], [1.95, 1.3, 0.0])
sq4_t = newEndPointsLineCurve("sq4_t", [1.95, 1.3, 0.0], [2.55, 1.3, 0.0])
sq4_r = newEndPointsLineCurve("sq4_r", [2.55, 1.3, 0.0], [2.55, 0.7, 0.0])

addCurveToInnerBoundary!(MY_project, sq4_b, "squareLoop4")
addCurveToInnerBoundary!(MY_project, sq4_l, "squareLoop4")
addCurveToInnerBoundary!(MY_project, sq4_t, "squareLoop4")
addCurveToInnerBoundary!(MY_project, sq4_r, "squareLoop4")

@info "All four inner square boundaries defined."

# ============================================================
#  REFINEMENT REGIONS
#  Centres and radii halved to match the new domain scale.
# ============================================================
#=
ref1 = newRefinementCenter("ref_sq1",  "smooth", [0.75, 3.5,  0.0], 0.15, 0.4)
ref2 = newRefinementCenter("ref_sq2",  "smooth", [2.25, 3.5,  0.0], 0.15, 0.4)
ref3 = newRefinementCenter("ref_sq3",  "smooth", [0.75, 1.0,  0.0], 0.15, 0.4)
ref4 = newRefinementCenter("ref_sq4",  "smooth", [2.25, 1.0,  0.0], 0.15, 0.4)
ref5 = newRefinementCenter("ref_vent", "smooth", [1.5,  4.5,  0.0], 0.15, 0.3)

addRefinementRegion!(MY_project, ref1)
addRefinementRegion!(MY_project, ref2)
addRefinementRegion!(MY_project, ref3)
addRefinementRegion!(MY_project, ref4)
addRefinementRegion!(MY_project, ref5)

@info "Refinement regions added."
=#
# ============================================================
#  VISUALISE & GENERATE
# ============================================================
plotProject!(MY_project, MODEL + GRID)

generate_mesh(MY_project)

@info "Mesh generation complete."

# ============================================================
#  TRIXI BOUNDARY CONDITIONS (reference — paste into sim file)
#  Boundary curve names are identical to the original mesh.
#
#  boundary_conditions = Dict(
#      :outer1   => boundary_condition_slip_wall,  # bottom
#      :outer2   => boundary_condition_slip_wall,  # right lower
#      :outer3   => boundary_condition_slip_wall,  # right notch h
#      :outer4   => boundary_condition_slip_wall,  # right notch v
#      :outer5   => boundary_condition_slip_wall,  # right notch h
#      :outer6   => boundary_condition_slip_wall,  # right upper
#      :outer7   => boundary_condition_slip_wall,  # top right
#      :vent     => boundary_condition_slip_wall,  # ← replace with inflow/outflow BC
#      :outer9   => boundary_condition_slip_wall,  # top left
#      :outer10  => boundary_condition_slip_wall,  # left upper
#      :outer11  => boundary_condition_slip_wall,  # left notch h
#      :outer12  => boundary_condition_slip_wall,  # left notch v
#      :outer13  => boundary_condition_slip_wall,  # left notch h
#      :outer14  => boundary_condition_slip_wall,  # left lower
#      :sq1_b    => boundary_condition_slip_wall,
#      :sq1_l    => boundary_condition_slip_wall,
#      :sq1_t    => boundary_condition_slip_wall,
#      :sq1_r    => boundary_condition_slip_wall,
#      :sq2_b    => boundary_condition_slip_wall,
#      :sq2_l    => boundary_condition_slip_wall,
#      :sq2_t    => boundary_condition_slip_wall,
#      :sq2_r    => boundary_condition_slip_wall,
#      :sq3_b    => boundary_condition_slip_wall,
#      :sq3_l    => boundary_condition_slip_wall,
#      :sq3_t    => boundary_condition_slip_wall,
#      :sq3_r    => boundary_condition_slip_wall,
#      :sq4_b    => boundary_condition_slip_wall,
#      :sq4_l    => boundary_condition_slip_wall,
#      :sq4_t    => boundary_condition_slip_wall,
#      :sq4_r    => boundary_condition_slip_wall,
#  )
#
#  REMEMBER: Update Blastwave_center and radius_initial in
#  your simulation file to match the new domain scale.
#  Original center (3, 4.5) → new center (1.5, 2.25)
#  Original radius 0.21875  → new radius 0.109375
# ============================================================