using VoronoiMeshes
using DelaunayTriangulation
using TensorsLite
using WriteVTK
using GLMakie

include("mesh_tools.jl")
using .MeshTools

const X_PERIOD = 1.0
const Y_PERIOD = 1.0

const MESH_PATTERN = r"^mesh_periodic_refined_nc(\d+)_vor\.vtu$"

# Independent random-point Lloyd relaxation needs a time budget that grows
# with nc, not a flat cap — a fixed budget either wastes time on small
# levels or, worse, silently cuts off large ones before they actually reach
# rtol. Scaled linearly (10 min per 1000 cells), with a floor for small
# levels and a 6-hour ceiling so the largest levels stay bounded; this is a
# simple calibration, not a rigorously derived formula.
const MINUTES_PER_1000_CELLS = 10.0
const MIN_MAX_TIME = 4.0     # package default, floor for tiny levels
const MAX_MAX_TIME = 360.0   # 6 hours, ceiling for huge levels

max_time_for(nc) = clamp(MINUTES_PER_1000_CELLS * nc / 1000, MIN_MAX_TIME, MAX_MAX_TIME)

# Every level is generated independently from random initial points (not a
# hex grid, not warm-started from another level) — warm-starting via cells +
# edge midpoints was tried and rejected: even a single bisection step
# reproduces strong hex-like regularity inherited from the parent mesh,
# regardless of how random that parent's own seed was.
function build_independent(nc)
    return VoronoiMesh(nc, X_PERIOD, Y_PERIOD, rtol=1e-5, max_iter=1_000_000, max_time=max_time_for(nc))
end

function build_level(outdir, mesh, i)
    nc = length(mesh.cells.position)
    label = "mesh_periodic_refined_nc$(nc)"
    return MeshTools.save_mesh_level(outdir, mesh, label, "p$i — $nc cells")
end

function main(base_cells, num_scales)
    outdir = MeshTools.run_outdir("refined_base$(base_cells)_n$(num_scales)")
    MeshTools.save_run_info(outdir, [
        "Script: build_set_refined_meshes_vtu.jl",
        "",
        "CLI parameters:",
        "  base_cells = $base_cells",
        "  num_scales = $num_scales",
        "",
        "Domain: periodic $(X_PERIOD) x $(Y_PERIOD)",
        "Density: uniform (no density function; VoronoiMesh default)",
        "Every level is independently generated from random initial points",
        "(rtol=1e-5, max_iter=1_000_000, max_time=$MINUTES_PER_1000_CELLS min per 1000 cells,",
        "floor $MIN_MAX_TIME min) -- not derived from one another and not warm-started,",
        "to avoid regularity inherited from bisection.",
        "Cell count quadruples each level (base_cells * 4^i).",
    ])

    rows = []
    for i in 0:(num_scales - 1)
        nc = base_cells * 4^i
        println("Scale p$i: $nc cells (independent random initial points, max_time=$(round(max_time_for(nc), digits=1)) min)...")
        mesh = build_independent(nc)
        push!(rows, build_level(outdir, mesh, i))
    end

    MeshTools.finalize_mesh_set(outdir, "refined", rows, MESH_PATTERN, MeshTools.numeric_sort_key(MESH_PATTERN))
    return nothing
end

const USAGE = """
Usage: julia --project=. build_set_refined_meshes_vtu.jl [base_cells] [num_scales]

Builds a series of independently-generated centroidal Voronoi meshes with
cell counts growing as base_cells * 4^i (not derived from one another, not
warm-started) — same cell-count ladder as build_set_regular_meshes.jl /
circular_refined, but each level seeded from fresh random points. Lloyd's
max_time scales with cell count (10 min per 1000 cells, 4 min floor, 6 hour
ceiling), so large levels can take a long time — up to 6h at nc=65536.

Arguments (all optional, positional):
  base_cells  Cell count at scale 0 (default 64).
  num_scales  Number of scales to build (default 4).
"""
MeshTools.handle_help(ARGS, USAGE)

base_cells = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 64
num_scales = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 4

main(base_cells, num_scales)
