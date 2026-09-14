using VoronoiMeshes
using DelaunayTriangulation
using TensorsLite
using WriteVTK
using GLMakie

include("mesh_tools.jl")
using .MeshTools

const X_PERIOD = 1.0
const Y_PERIOD = 1.0

const MESH_PATTERN = r"^mesh_periodic_regular_nc(\d+)_vor\.vtu$"

function build_level(outdir, mesh)
    nc = length(mesh.cells.position)
    label = "mesh_periodic_regular_nc$(nc)"
    return MeshTools.save_mesh_level(outdir, mesh, label, "$nc cells")
end

function main(nc_ref, num_levels)
    outdir = MeshTools.run_outdir("regular_nc$(nc_ref)_L$(num_levels)")
    MeshTools.save_run_info(outdir, [
        "Script: build_set_regular_meshes.jl",
        "",
        "CLI parameters:",
        "  nc_ref = $nc_ref",
        "  num_levels = $num_levels",
        "",
        "Domain: periodic $(X_PERIOD) x $(Y_PERIOD)",
        "Density: uniform (no density function; VoronoiMesh default)",
        "Each level refines from the previous one's cells + edge midpoints as new generators.",
    ])

    mesh, dc = MeshTools.build_hex_reference(nc_ref, X_PERIOD, Y_PERIOD)
    println("Level 0: creating regular hex mesh (reference nc≈$nc_ref, dc=$(round(dc, digits=4)))...")
    rows = [build_level(outdir, mesh)]

    for i in 1:num_levels
        nc_prev = length(mesh.cells.position)
        println("Level $i: refining from $nc_prev cells (cells + edge midpoints as generators)...")
        generators = vcat(mesh.cells.position, mesh.edges.position)
        mesh = VoronoiMesh(generators, X_PERIOD, Y_PERIOD)
        push!(rows, build_level(outdir, mesh))
    end

    MeshTools.finalize_mesh_set(outdir, "regular", rows, MESH_PATTERN, MeshTools.numeric_sort_key(MESH_PATTERN))
    return nothing
end

const USAGE = """
Usage: julia --project=. build_set_regular_meshes.jl [nc_ref] [num_levels]

Builds a series of regular Voronoi meshes by successive refinement, starting
from a regular hex mesh and quadrupling the cell count each level.

Arguments (all optional, positional):
  nc_ref      Target cell count for the Level 0 mesh (default 16).
  num_levels  Number of refinement levels after Level 0 (default 4).
"""
MeshTools.handle_help(ARGS, USAGE)

nc_ref     = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 16
num_levels = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 4

main(nc_ref, num_levels)
