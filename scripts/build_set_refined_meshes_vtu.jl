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
    outdir = MeshTools.run_outdir("refined_base$(base_cells)_n$(num_scales)_i$(ini_scale)")
    MeshTools.save_run_info(outdir, [
        "Script: build_set_refined_meshes_vtu.jl",
        "",
        "CLI parameters:",
        "  base_cells = $base_cells",
        "  num_scales = $num_scales",
        "  ini_scale = $ini_scale",
        "",
        "Domain: periodic $(X_PERIOD) x $(Y_PERIOD)",
        "Density: uniform (no density function; VoronoiMesh default)",
        "Each scale is an independently-generated centroidal mesh (not derived from the others);",
        "rtol=1e-5, max_iter=100000, max_time=10.0 min",
    ])

    rows = []
    for i in ini_scale:(ini_scale+num_scales-1)
        num_cells = base_cells * (4^i)
        println("Scale p$i: creating centroidal mesh ($num_cells cells)...")

        # max_time raised from the 4-minute default: an early time-cap cutoff
        # leaves a still-shifting generator set that can crash mesh construction.
        mesh = VoronoiMesh(num_cells, X_PERIOD, Y_PERIOD, rtol=1e-5, max_iter=100000, max_time=10.0)

        label = "mesh_periodic_refined_nc$(num_cells)"
        push!(rows, MeshTools.save_mesh_level(outdir, mesh, label, "p$i — $num_cells cells"))
    end

    MeshTools.finalize_mesh_set(outdir, "refined", rows, MESH_PATTERN, MeshTools.numeric_sort_key(MESH_PATTERN))
    return nothing
end

const USAGE = """
Usage: julia --project=. build_set_refined_meshes_vtu.jl [base_cells] [num_scales] [ini_scale]

Builds a series of independently-generated centroidal Voronoi meshes with
cell counts growing as base_cells * 4^i (not derived from one another) —
same cell-count ladder as build_set_regular_meshes.jl / circular_refined.

Arguments (all optional, positional):
  base_cells  Cell count at ini_scale (default 64).
  num_scales  Number of scales to build (default 4).
  ini_scale   Starting power of 4 (default 0). Cell counts:
              base_cells*4^ini_scale, ..., base_cells*4^(ini_scale+num_scales-1)
              — defaults give 64, 256, 1024, 4096.
"""
MeshTools.handle_help(ARGS, USAGE)

base_cells = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 64
num_scales = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 4
ini_scale = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 0

main(base_cells, num_scales, ini_scale)
