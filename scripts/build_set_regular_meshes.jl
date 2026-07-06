# build_set_regular_meshes.jl
#
# Builds a series of regular Voronoi meshes by successive refinement starting
# from a regular hexagonal mesh.
#
#   Level 0 : regular hex mesh from create_planar_hex_mesh (~nc_ref cells)
#   Level i : uses cell centres + edge midpoints of level i-1 as generators,
#             then runs Lloyd relaxation to convergence.
#             This roughly quadruples the cell count at each level.
#
# Usage:
#   julia --project=. build_set_regular_meshes.jl [nc_ref] [num_levels]
#
# Defaults: nc_ref=16, num_levels=4

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
    outdir = "output"
    mkpath(outdir)

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

nc_ref     = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 16
num_levels = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 4

main(nc_ref, num_levels)
