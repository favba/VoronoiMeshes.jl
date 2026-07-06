# build_set_refined_meshes_vtu.jl
#
# Builds a series of centroidal Voronoi meshes with increasing cell counts,
# each generated independently from random initial points and converged via
# Lloyd's algorithm. Cell counts grow as base_cells * 2^i.
#
# Usage:
#   julia --project=. build_set_refined_meshes_vtu.jl [base_cells] [num_scales] [ini_scale]
#
# Defaults: base_cells=16, num_scales=11, ini_scale=0
#
# Cell counts: base_cells*2^ini_scale, base_cells*2^(ini_scale+1), ...,
# base_cells*2^(ini_scale+num_scales-1). E.g. defaults give 16, 32, 64, ..., 16384.
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

const MESH_PATTERN = r"^mesh_periodic_refined_nc(\d+)_vor\.vtu$"

function main(base_cells, num_scales, ini_scale)
    outdir = "output"
    mkpath(outdir)

    rows = []
    for i in ini_scale:(ini_scale+num_scales-1)
        num_cells = base_cells * (2^i)
        println("Scale p$i: creating centroidal mesh ($num_cells cells)...")

        # Independently-generated centroidal Voronoi mesh (converged via Lloyd's
        # algorithm), not derived from the previous scale's generators.
        mesh = VoronoiMesh(num_cells, X_PERIOD, Y_PERIOD, rtol=1e-5, max_iter=100000)

        label = "mesh_periodic_refined_nc$(num_cells)"
        push!(rows, MeshTools.save_mesh_level(outdir, mesh, label, "p$i — $num_cells cells"))
    end

    MeshTools.finalize_mesh_set(outdir, "refined", rows, MESH_PATTERN, MeshTools.numeric_sort_key(MESH_PATTERN))
    return nothing
end

base_cells = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 16
num_scales = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 11
ini_scale = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 0

main(base_cells, num_scales, ini_scale)
