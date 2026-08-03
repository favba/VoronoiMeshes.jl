# build_set_localref_meshes.jl
#
# Builds a series of locally-refined variable-resolution meshes: density is
# RIDGE_RATIO times higher along a diagonal band (LINE_SLOPE/LINE_INTERCEPT,
# same line as build_set_irregular_meshes.jl) than in the background,
# smoothly (Gaussian) tapering over RIDGE_SIGMA. Refinement quadruples nc
# each level, doubling both the ridge and background resolution while the
# fixed density-function shape keeps their ratio constant.
#
# RIDGE_RATIO = 16 comes from Du, Faber & Gunzburger (1999), "Centroidal
# Voronoi Tessellations: Applications and Algorithms" (SIAM Review): for a
# CVT with density ρ(x) in R^d, the asymptotic point density satisfies
# m(x) ∝ ρ(x)^(d/(d+2)); in 2D, m(x) ∝ ρ(x)^(1/2), and since cell area ∝
# 1/m(x) and dc ∝ sqrt(area), dc(x) ∝ ρ(x)^(-1/4). Solving R^(-1/4) = 1/2
# (target: ridge dc is half the background dc, i.e. 2x resolution) gives
# R = 2^4 = 16. This is an asymptotic result, not an exact guarantee at
# finite resolution, so each level's actual ridge/background dc ratio is
# measured directly (report_resolution_ratio) rather than trusted blindly.
#
# Usage:
#   julia --project=. build_set_localref_meshes.jl [nc_ref] [num_levels]
#
# Defaults: nc_ref=80, num_levels=5

using VoronoiMeshes
using DelaunayTriangulation
using TensorsLite
using WriteVTK
using GLMakie
using Statistics: mean

include("mesh_tools.jl")
using .MeshTools

const X_PERIOD = 1.0
const Y_PERIOD = 1.0

# Same diagonal line as build_set_irregular_meshes.jl.
const LINE_SLOPE = 0.3
const LINE_INTERCEPT = 0.3

const RIDGE_SIGMA = 0.06   # smooth-transition width (Gaussian sigma); tune after visual/measured check
const RIDGE_RATIO = 16.0   # see header comment for derivation

# Density-weighted Lloyd relaxation with this much contrast converges far
# more slowly than the uniform case: even at nc~50-200 it exhausts a
# 20_000-iteration budget without reaching rtol=1e-3, and the package's
# default max_time (4 minutes) cuts it off before that too. Loosen the
# tolerance and budget much more time/iterations than the package defaults —
# expect this script to run far longer than the other build scripts,
# especially at the finer levels.
const LOCALREF_RTOL = 1.0e-3
const LOCALREF_MAX_ITER = 50_000
const LOCALREF_MAX_TIME = 60.0   # minutes, per level

ridge_density(p) = 1 + (RIDGE_RATIO - 1) *
    exp(-0.5 * (MeshTools.line_distance(p, LINE_SLOPE, LINE_INTERCEPT) / RIDGE_SIGMA)^2)

const MESH_PATTERN = r"^mesh_periodic_localref_nc(\d+)_vor\.vtu$"

# Measures and prints the actual ridge-vs-background cell-diameter ratio —
# validates the intended resolution contrast regardless of how RIDGE_RATIO
# was chosen (see header comment: the density-ratio-to-dc-ratio relationship
# is asymptotic, not exact, so it's checked directly rather than trusted).
function report_resolution_ratio(mesh)
    d = MeshTools.cell_diameter(mesh)
    on_ridge = [MeshTools.line_distance(p, LINE_SLOPE, LINE_INTERCEPT) < RIDGE_SIGMA
                for p in mesh.cells.position]
    ratio = mean(d[on_ridge]) / mean(d[.!on_ridge])
    println("    ridge/background dc ratio = $(round(ratio, digits=3)) (target 0.5)")
end

function build_level(outdir, mesh)
    nc = length(mesh.cells.position)
    label = "mesh_periodic_localref_nc$(nc)"
    row = MeshTools.save_mesh_level(outdir, mesh, label, "$nc cells")
    report_resolution_ratio(mesh)
    return row
end

function main(nc_ref, num_levels)
    outdir = "output"
    mkpath(outdir)

    mesh, dc = MeshTools.build_hex_reference(
        nc_ref, X_PERIOD, Y_PERIOD;
        density=ridge_density, rtol=LOCALREF_RTOL, max_iter=LOCALREF_MAX_ITER, max_time=LOCALREF_MAX_TIME,
    )
    println("Level 0: creating hex mesh under ridge density (reference nc≈$nc_ref, dc=$(round(dc, digits=4)))...")
    rows = [build_level(outdir, mesh)]

    for i in 1:(num_levels - 1)
        nc_prev = length(mesh.cells.position)
        println("Level $i: refining from $nc_prev cells (cells + edge midpoints as generators)...")
        generators = vcat(mesh.cells.position, mesh.edges.position)
        mesh = VoronoiMesh(
            generators, X_PERIOD, Y_PERIOD;
            density=ridge_density, rtol=LOCALREF_RTOL, max_iter=LOCALREF_MAX_ITER, max_time=LOCALREF_MAX_TIME,
        )
        push!(rows, build_level(outdir, mesh))
    end

    MeshTools.finalize_mesh_set(outdir, "localref", rows, MESH_PATTERN, MeshTools.numeric_sort_key(MESH_PATTERN))
    return nothing
end

nc_ref     = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 80
num_levels = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 5

main(nc_ref, num_levels)
