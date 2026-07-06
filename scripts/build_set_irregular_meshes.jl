# build_set_irregular_meshes.jl
#
# Generates a series of nc-cell Voronoi meshes with increasing irregularity.
#   Level 0 : fully converged centroidal mesh (most regular)
#   Levels 1..num_levels : each level takes the previous level's generator
#              points, adds growing random noise, then runs PERTURB_ITERS
#              Lloyd steps — enough to fix degenerate cells without erasing
#              the introduced irregularity.
#
# Usage:
#   julia --project=. build_set_irregular_meshes.jl [nc] [num_levels] [base_strength]
#   julia --project=. build_set_irregular_meshes.jl [nc] 0 [strength]   # single perturbed mesh
#
# Defaults: nc=80, num_levels=5, base_strength=0.15
#
# base_strength is the perturbation amplitude at level 1, expressed as a
# fraction of the average cell spacing dc ≈ 1/√nc.  It scales linearly:
# level i uses strength = base_strength * i.
#
# Use num_levels=0 to build only the Level 0 mesh plus a single perturbation
# at exactly base_strength (no multi-level sweep).
#
# Note: PNG export uses GLMakie. On headless servers swap to CairoMakie
# (add it to the environment and replace `using GLMakie` below).

using VoronoiMeshes
using DelaunayTriangulation
using TensorsLite
using WriteVTK
using GLMakie

include("mesh_tools.jl")
using .MeshTools

const X_PERIOD = 1.0
const Y_PERIOD = 1.0

const PERTURB_ITERS = 5

# Perturbation band: only cells within BAND_WIDTH of the line y = LINE_SLOPE*x + LINE_INTERCEPT
# are perturbed. Cells outside the band keep their positions unchanged.
const LINE_SLOPE = 0.3
const LINE_INTERCEPT = 0.3
const BAND_WIDTH = 0.12

# Per-cell displacement is clamped to this fraction of dc, in absolute terms —
# not a multiple of strength*dc, since that would keep scaling up with strength
# and never actually bound anything. Without this cap, a displacement
# comparable to a cell's own radius (not just a rare tail draw) is what
# actually causes overlap/inversion (obtuse Delaunay triangles / localized
# defects), and that starts happening for *typical* draws once strength alone
# gets large — clamping the Gaussian's sigma multiplier doesn't help since it
# scales right along with strength.
const MAX_DISPLACEMENT_FRAC = 0.3

# If a perturbation still produces obtuse Delaunay triangles (circumcenter
# outside triangle) despite the clamp above, redraw the whole band's random
# perturbation from scratch (fresh random numbers) rather than trying to
# Lloyd-relax the defect away — up to this many attempts.
const MAX_REDRAWS = 10

const MESH_PATTERN = r"^mesh_periodic_irregular_nc\d+_d([\d.]+)_vor\.vtu$"

# Obtuse-triangle count, printed after the shared metrics summary (already
# printed by MeshTools.save_mesh_level) — specific to this script's
# redraw-on-defect loop.
function print_obtuse_triangles(mesh)
    n_obtuse = length(find_obtuse_triangles(mesh))
    n_tri = length(mesh.vertices.cells)
    println("    obtuse triangles = $n_obtuse / $n_tri")
end

function perturb_positions(mesh, strength)
    pos = mesh.cells.position
    nc = length(pos)
    dc = sqrt(X_PERIOD * Y_PERIOD / nc)
    max_disp = MAX_DISPLACEMENT_FRAC * dc

    # Perpendicular distance from (x,y) to line y = LINE_SLOPE*x + LINE_INTERCEPT,
    # i.e. -LINE_SLOPE*x + y - LINE_INTERCEPT = 0
    norm_factor = sqrt(1 + LINE_SLOPE^2)

    x_new = copy(pos.x)
    y_new = copy(pos.y)
    for c in 1:nc
        p = pos[c]
        dist = abs(p.y - LINE_SLOPE * p.x - LINE_INTERCEPT) / norm_factor
        if dist < BAND_WIDTH
            dx = clamp(strength * dc * randn(), -max_disp, max_disp)
            dy = clamp(strength * dc * randn(), -max_disp, max_disp)
            x_new[c] = mod(p.x + dx, X_PERIOD)
            y_new[c] = mod(p.y + dy, Y_PERIOD)
        end
    end
    return VecArray(x=x_new, y=y_new)
end

function build_reference_mesh(nc, outdir)
    mesh, dc = MeshTools.build_hex_reference(nc, X_PERIOD, Y_PERIOD)
    actual_nc = length(mesh.cells.position)
    println("Level 0: creating regular hex mesh (reference nc≈$nc, actual=$actual_nc, dc=$(round(dc, digits=4)))...")

    label = "mesh_periodic_irregular_nc$(actual_nc)_d0.0"
    row = MeshTools.save_mesh_level(outdir, mesh, label)
    print_obtuse_triangles(mesh)

    return mesh, actual_nc, row
end

function perturb_level(prev_mesh, strength, level, actual_nc, outdir)
    println("Level $level: perturbing (strength=$(round(strength, digits=3))), $PERTURB_ITERS Lloyd iterations...")

    mesh = nothing
    n_obtuse = 0
    redraws = 0
    while true
        generators = perturb_positions(prev_mesh, strength)
        mesh = VoronoiMesh(generators, X_PERIOD, Y_PERIOD; max_iter=PERTURB_ITERS)
        n_obtuse = length(find_obtuse_triangles(mesh))
        (n_obtuse == 0 || redraws >= MAX_REDRAWS) && break
        redraws += 1
    end
    if redraws > 0
        status = n_obtuse == 0 ? "resolved" : "still $n_obtuse remaining after $MAX_REDRAWS attempts"
        println("  Redrew perturbation $redraws time(s) to avoid obtuse triangles ($status).")
    end

    label = "mesh_periodic_irregular_nc$(actual_nc)_d$(round(strength, digits=3))"
    row = MeshTools.save_mesh_level(outdir, mesh, label)
    print_obtuse_triangles(mesh)

    return mesh, row
end

function main(nc, num_levels, base_strength)
    outdir = "output"
    mkpath(outdir)

    mesh, actual_nc, row0 = build_reference_mesh(nc, outdir)
    rows = [row0]

    if num_levels == 0
        # Single mesh, single perturbation at exactly base_strength.
        _, row = perturb_level(mesh, base_strength, 1, actual_nc, outdir)
        push!(rows, row)
    else
        prev_mesh = mesh
        for i in 1:num_levels
            strength = base_strength * i
            prev_mesh, row = perturb_level(prev_mesh, strength, i, actual_nc, outdir)
            push!(rows, row)
        end
    end

    MeshTools.finalize_mesh_set(outdir, "irregular", rows, MESH_PATTERN, MeshTools.numeric_sort_key(MESH_PATTERN))
    return nothing
end

nc = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 80
num_levels = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 5
base_strength = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 0.15
base_strength < 0 && error("base_strength must be >= 0 (got $base_strength)")

main(nc, num_levels, base_strength)
