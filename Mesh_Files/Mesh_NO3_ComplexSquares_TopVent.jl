using HOHQMesh, GLMakie

# ============================================================
#  Project setup
# ============================================================
MY_project = newProject("Mesh_NO3_ComplexSquares_TopVent", "Mesh_Files/Meshes_Folder")

setPolynomialOrder!(MY_project, 3)
setMeshFileFormat!(MY_project, "ABAQUS")

addBackgroundGrid!(MY_project, [0.25, 0.25, 0.0])

# ============================================================
#  OUTER BOUNDARY
#
#  Domain is a 6 × 9 rectangle (x: 0→6, y: 0→9) with:
#    - A notch cut into the LEFT  wall at y = 3.0 → 4.5
#    - A notch cut into the RIGHT wall at y = 3.0 → 4.5
#    - A vent OPENING at the TOP  from x = 2.5 → 3.5
#      (open boundary — no slip wall on that segment)
#
#  Traversal order: counter-clockwise starting at origin.
# ============================================================

# --- Bottom edge ---
outer1  = newEndPointsLineCurve("outer1",  [0.0, 0.0, 0.0], [6.0, 0.0, 0.0])

# --- Right wall: bottom section ---
outer2  = newEndPointsLineCurve("outer2",  [6.0, 0.0, 0.0], [6.0, 3.0, 0.0])

# --- Right notch (indents inward, i.e. toward x = 5.5) ---
outer3  = newEndPointsLineCurve("outer3",  [6.0, 3.0, 0.0], [5.5, 3.0, 0.0])
outer4  = newEndPointsLineCurve("outer4",  [5.5, 3.0, 0.0], [5.5, 4.5, 0.0])
outer5  = newEndPointsLineCurve("outer5",  [5.5, 4.5, 0.0], [6.0, 4.5, 0.0])

# --- Right wall: top section ---
outer6  = newEndPointsLineCurve("outer6",  [6.0, 4.5, 0.0], [6.0, 9.0, 0.0])

# --- Top edge: right portion (x = 6 → 3.5) ---
outer7  = newEndPointsLineCurve("outer7",  [6.0, 9.0, 0.0], [3.5, 9.0, 0.0])

# --- Vent opening at top (x = 3.5 → 2.5): open boundary ---
#     Give it a distinct name so it can receive a different
#     boundary condition (e.g. inflow / outflow) in Trixi.
outer8  = newEndPointsLineCurve("vent",    [3.5, 9.0, 0.0], [2.5, 9.0, 0.0])

# --- Top edge: left portion (x = 2.5 → 0) ---
outer9  = newEndPointsLineCurve("outer9",  [2.5, 9.0, 0.0], [0.0, 9.0, 0.0])

# --- Left wall: top section ---
outer10 = newEndPointsLineCurve("outer10", [0.0, 9.0, 0.0], [0.0, 4.5, 0.0])

# --- Left notch (indents inward, i.e. toward x = 0.5) ---
outer11 = newEndPointsLineCurve("outer11", [0.0, 4.5, 0.0], [0.5, 4.5, 0.0])
outer12 = newEndPointsLineCurve("outer12", [0.5, 4.5, 0.0], [0.5, 3.0, 0.0])
outer13 = newEndPointsLineCurve("outer13", [0.5, 3.0, 0.0], [0.0, 3.0, 0.0])

# --- Left wall: bottom section ---
outer14 = newEndPointsLineCurve("outer14", [0.0, 3.0, 0.0], [0.0, 0.0, 0.0])

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
#  Square centres and half-widths:
#    sq1: centre (1.5, 7.0), half = 0.6  → x: 0.9→2.1, y: 6.4→7.6
#    sq2: centre (4.5, 7.0), half = 0.6  → x: 3.9→5.1, y: 6.4→7.6
#    sq3: centre (1.5, 2.0), half = 0.6  → x: 0.9→2.1, y: 1.4→2.6
#    sq4: centre (4.5, 2.0), half = 0.6  → x: 3.9→5.1, y: 1.4→2.6
# ============================================================

# --- Square hole 1 (top-left) ---
sq1_b = newEndPointsLineCurve("sq1_b", [2.1, 6.4, 0.0], [0.9, 6.4, 0.0])  # bottom, right→left
sq1_l = newEndPointsLineCurve("sq1_l", [0.9, 6.4, 0.0], [0.9, 7.6, 0.0])  # left,   bottom→top
sq1_t = newEndPointsLineCurve("sq1_t", [0.9, 7.6, 0.0], [2.1, 7.6, 0.0])  # top,    left→right
sq1_r = newEndPointsLineCurve("sq1_r", [2.1, 7.6, 0.0], [2.1, 6.4, 0.0])  # right,  top→bottom

addCurveToInnerBoundary!(MY_project, sq1_b, "squareLoop1")
addCurveToInnerBoundary!(MY_project, sq1_l, "squareLoop1")
addCurveToInnerBoundary!(MY_project, sq1_t, "squareLoop1")
addCurveToInnerBoundary!(MY_project, sq1_r, "squareLoop1")

# --- Square hole 2 (top-right) ---
sq2_b = newEndPointsLineCurve("sq2_b", [5.1, 6.4, 0.0], [3.9, 6.4, 0.0])
sq2_l = newEndPointsLineCurve("sq2_l", [3.9, 6.4, 0.0], [3.9, 7.6, 0.0])
sq2_t = newEndPointsLineCurve("sq2_t", [3.9, 7.6, 0.0], [5.1, 7.6, 0.0])
sq2_r = newEndPointsLineCurve("sq2_r", [5.1, 7.6, 0.0], [5.1, 6.4, 0.0])

addCurveToInnerBoundary!(MY_project, sq2_b, "squareLoop2")
addCurveToInnerBoundary!(MY_project, sq2_l, "squareLoop2")
addCurveToInnerBoundary!(MY_project, sq2_t, "squareLoop2")
addCurveToInnerBoundary!(MY_project, sq2_r, "squareLoop2")

# --- Square hole 3 (bottom-left) ---
sq3_b = newEndPointsLineCurve("sq3_b", [2.1, 1.4, 0.0], [0.9, 1.4, 0.0])
sq3_l = newEndPointsLineCurve("sq3_l", [0.9, 1.4, 0.0], [0.9, 2.6, 0.0])
sq3_t = newEndPointsLineCurve("sq3_t", [0.9, 2.6, 0.0], [2.1, 2.6, 0.0])
sq3_r = newEndPointsLineCurve("sq3_r", [2.1, 2.6, 0.0], [2.1, 1.4, 0.0])

addCurveToInnerBoundary!(MY_project, sq3_b, "squareLoop3")
addCurveToInnerBoundary!(MY_project, sq3_l, "squareLoop3")
addCurveToInnerBoundary!(MY_project, sq3_t, "squareLoop3")
addCurveToInnerBoundary!(MY_project, sq3_r, "squareLoop3")

# --- Square hole 4 (bottom-right) ---
sq4_b = newEndPointsLineCurve("sq4_b", [5.1, 1.4, 0.0], [3.9, 1.4, 0.0])
sq4_l = newEndPointsLineCurve("sq4_l", [3.9, 1.4, 0.0], [3.9, 2.6, 0.0])
sq4_t = newEndPointsLineCurve("sq4_t", [3.9, 2.6, 0.0], [5.1, 2.6, 0.0])
sq4_r = newEndPointsLineCurve("sq4_r", [5.1, 2.6, 0.0], [5.1, 1.4, 0.0])

addCurveToInnerBoundary!(MY_project, sq4_b, "squareLoop4")
addCurveToInnerBoundary!(MY_project, sq4_l, "squareLoop4")
addCurveToInnerBoundary!(MY_project, sq4_t, "squareLoop4")
addCurveToInnerBoundary!(MY_project, sq4_r, "squareLoop4")

@info "All four inner square boundaries defined."

# ============================================================
#  REFINEMENT REGIONS
#  Place smoothing centres near each obstacle and the vent
#  to get better resolution around features of interest.
# ============================================================

ref1 = newRefinementCenter("ref_sq1",  "smooth", [1.5, 7.0, 0.0], 0.3, 0.8)
ref2 = newRefinementCenter("ref_sq2",  "smooth", [4.5, 7.0, 0.0], 0.3, 0.8)
ref3 = newRefinementCenter("ref_sq3",  "smooth", [1.5, 2.0, 0.0], 0.3, 0.8)
ref4 = newRefinementCenter("ref_sq4",  "smooth", [4.5, 2.0, 0.0], 0.3, 0.8)
ref5 = newRefinementCenter("ref_vent", "smooth", [3.0, 9.0, 0.0], 0.3, 0.6)

addRefinementRegion!(MY_project, ref1)
addRefinementRegion!(MY_project, ref2)
addRefinementRegion!(MY_project, ref3)
addRefinementRegion!(MY_project, ref4)
addRefinementRegion!(MY_project, ref5)

@info "Refinement regions added."

# ============================================================
#  VISUALISE & GENERATE
# ============================================================
plotProject!(MY_project, MODEL + GRID)

generate_mesh(MY_project)

@info "Mesh generation complete."

# ============================================================
#  TRIXI BOUNDARY CONDITIONS (reference — paste into sim file)
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
# ============================================================