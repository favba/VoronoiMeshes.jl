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

# create_planar_hex_mesh doesn't always land exactly on an arbitrary target
# nc (it rounds to the nearest valid hex row/column count), so Level 0 is
# bootstrapped from a small exact mesh and bisected up — same approach as
# build_set_circular_refined_meshes.jl — guaranteeing nc_ref lands exactly
# whenever it's reachable by quadrupling from BOOTSTRAP_NC.
const BOOTSTRAP_NC = 16

function refine_step(mesh)
    generators = vcat(mesh.cells.position, mesh.edges.position)
    return VoronoiMesh(generators, X_PERIOD, Y_PERIOD)
end

function build_level0(nc_ref)
    mesh, _ = MeshTools.build_hex_reference(min(nc_ref, BOOTSTRAP_NC), X_PERIOD, Y_PERIOD)
    while length(mesh.cells.position) < nc_ref
        mesh = refine_step(mesh)
    end
    return mesh
end

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

    mesh = build_level0(nc_ref)
    println("Level 0: regular hex mesh (reference nc≈$nc_ref, actual=$(length(mesh.cells.position)))...")
    rows = [build_level(outdir, mesh)]

    for i in 1:(num_levels - 1)
        nc_prev = length(mesh.cells.position)
        println("Level $i: refining from $nc_prev cells (cells + edge midpoints as generators)...")
        mesh = refine_step(mesh)
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
  nc_ref      Target cell count for Level 0 (default 64).
  num_levels  Number of levels, Level 0 + refinements (default 4).
"""
MeshTools.handle_help(ARGS, USAGE)

nc_ref     = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 64
num_levels = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 4

main(nc_ref, num_levels)
