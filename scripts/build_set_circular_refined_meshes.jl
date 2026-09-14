using VoronoiMeshes
using DelaunayTriangulation
using TensorsLite
using WriteVTK
using GLMakie

include("mesh_tools.jl")
using .MeshTools

const X_PERIOD = 1.0
const Y_PERIOD = 1.0

# Off-center refinement target; periodicity holds as long as
# buffer_radius < min(lx,ly)/2, regardless of where CENTER sits.
const CENTER = (2/3) * X_PERIOD * 𝐢 + (2/3) * Y_PERIOD * 𝐣
const INNER_RADIUS = 1 / 6
const TRANSITION_WIDTH = INNER_RADIUS
const BUFFER_RADIUS = INNER_RADIUS + TRANSITION_WIDTH

# dc(x) ∝ ρ(x)^(-1/4) in 2D (Du, Faber & Gunzburger 1999), so R=16 targets a
# 2x inner/outer resolution ratio.
const DENSITY_RATIO = 16.0

# Density-weighted Lloyd at this contrast converges far more slowly than
# uniform relaxation; loosened well beyond the package defaults.
const CIRCULAR_RTOL = 1.0e-3
const CIRCULAR_MAX_ITER = 50_000
const CIRCULAR_MAX_TIME = 60.0   # minutes, per level

const circular_density = VoronoiMeshes.circular_refinement_function(
    X_PERIOD, Y_PERIOD;
    center = CENTER, inner_radius = INNER_RADIUS, buffer_radius = BUFFER_RADIUS,
    inner_density = DENSITY_RATIO, outer_density = 1.0,
)

const MESH_PATTERN = r"^mesh_periodic_circular_refined_nc(\d+)_vor\.vtu$"

# Converging Level 0 directly at nc_ref is the slowest step (the full density
# contrast has to develop from scratch), so bootstrap from a small mesh
# already converged under the same density and bisect it up instead.
const BOOTSTRAP_NC = 16

function refine_step(mesh)
    generators = vcat(mesh.cells.position, mesh.edges.position)
    return VoronoiMesh(
        generators, X_PERIOD, Y_PERIOD;
        density = circular_density, rtol = CIRCULAR_RTOL, max_iter = CIRCULAR_MAX_ITER, max_time = CIRCULAR_MAX_TIME,
    )
end

function build_level0(nc_ref)
    mesh, _ = MeshTools.build_hex_reference(
        min(nc_ref, BOOTSTRAP_NC), X_PERIOD, Y_PERIOD;
        density = circular_density, rtol = CIRCULAR_RTOL, max_iter = CIRCULAR_MAX_ITER, max_time = CIRCULAR_MAX_TIME,
    )
    while length(mesh.cells.position) < nc_ref
        nc_prev = length(mesh.cells.position)
        println("Bootstrap: refining from $nc_prev cells toward nc_ref≈$nc_ref...")
        mesh = refine_step(mesh)
    end
    return mesh
end

function build_level(outdir, mesh)
    nc = length(mesh.cells.position)
    label = "mesh_periodic_circular_refined_nc$(nc)"
    return MeshTools.save_mesh_level(outdir, mesh, label, "$nc cells"; region = (CENTER, INNER_RADIUS, BUFFER_RADIUS, DENSITY_RATIO))
end

function main(nc_ref, num_levels)
    outdir = MeshTools.run_outdir("circular_refined_nc$(nc_ref)_L$(num_levels)")
    MeshTools.save_run_info(outdir, [
        "Script: build_set_circular_refined_meshes.jl",
        "",
        "CLI parameters:",
        "  nc_ref = $nc_ref",
        "  num_levels = $num_levels",
        "",
        "Domain: periodic $(X_PERIOD) x $(Y_PERIOD)",
        "",
        "Density function: VoronoiMeshes.circular_refinement_function",
        "  center = ($(CENTER.x), $(CENTER.y))",
        "  inner_radius = $INNER_RADIUS",
        "  buffer_radius = $BUFFER_RADIUS  (transition_width = $TRANSITION_WIDTH)",
        "  inner_density = $DENSITY_RATIO",
        "  outer_density = 1.0",
        "",
        "Lloyd relaxation: rtol=$CIRCULAR_RTOL, max_iter=$CIRCULAR_MAX_ITER, max_time=$CIRCULAR_MAX_TIME min/level",
        "Level 0 bootstrapped from BOOTSTRAP_NC=$BOOTSTRAP_NC, bisected up to nc_ref.",
    ])

    mesh = build_level0(nc_ref)
    println("Level 0: circular-density mesh (reference nc≈$nc_ref, actual=$(length(mesh.cells.position)))...")
    rows = [build_level(outdir, mesh)]

    for i in 1:(num_levels - 1)
        nc_prev = length(mesh.cells.position)
        println("Level $i: refining from $nc_prev cells (cells + edge midpoints as generators)...")
        mesh = refine_step(mesh)
        push!(rows, build_level(outdir, mesh))
    end

    MeshTools.finalize_mesh_set(outdir, "circular_refined", rows, MESH_PATTERN, MeshTools.numeric_sort_key(MESH_PATTERN))
    return nothing
end

const USAGE = """
Usage: julia --project=. build_set_circular_refined_meshes.jl [nc_ref] [num_levels]

Builds a series of locally-refined variable-resolution meshes around an
off-center circular refinement region, refining (quadrupling nc) each level
at fixed transition sharpness and density ratio.

Arguments (all optional, positional):
  nc_ref      Target cell count for the Level 0 mesh (default 64).
  num_levels  Number of levels, Level 0 + refinements (default 4).

Example (DG-test scale, finest level ~30k cells):
  julia --project=. build_set_circular_refined_meshes.jl 32 6
"""
MeshTools.handle_help(ARGS, USAGE)

nc_ref     = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 64
num_levels = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 4

main(nc_ref, num_levels)
