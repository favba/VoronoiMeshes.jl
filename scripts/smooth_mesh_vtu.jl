using VoronoiMeshes
using DelaunayTriangulation
using TensorsLite
using WriteVTK
using ReadVTK
using GLMakie

include("mesh_tools.jl")
using .MeshTools

const USAGE = """
Usage: julia --project=. smooth_mesh_vtu.jl <mesh_vor.vtu> [max_time] [rtol]

Reads a previously-saved mesh and applies further Lloyd relaxation,
warm-started from its own cell positions (same cell count, no
refinement/upsampling) — for touching up a mesh that stopped early (e.g.
hit a time cap) without rebuilding it from scratch. Saves the smoothed
mesh (VTU) and per-metric plots to its own output/ subfolder.

Arguments:
  mesh_vor.vtu  Path to a "*_vor.vtu" file (its matching "*_tri.vtu" is
                found automatically alongside it).
  max_time      Lloyd relaxation time budget, in minutes (default 10.0).
  rtol          Lloyd relative tolerance (default 1e-5).
"""
MeshTools.handle_help(ARGS, USAGE)

isempty(ARGS) && error(USAGE)
vtu_path = ARGS[1]
max_time = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 10.0
rtol = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 1e-5

tri_path = replace(vtu_path, "_vor.vtu" => "_tri.vtu")
mesh = VoronoiMeshes.read_from_vtu(vtu_path, tri_path)
nc = length(mesh.cells.position)
xp, yp = mesh.x_period, mesh.y_period
println("Loaded $vtu_path: nc=$nc, period=($xp, $yp)")

before_obtuse, before_tri = MeshTools.obtuse_triangle_count(mesh)
println("Before smoothing: obtuse triangles = $before_obtuse / $before_tri")

println("Smoothing (rtol=$rtol, max_time=$max_time min)...")
smoothed = VoronoiMesh(mesh.cells.position, xp, yp; rtol=rtol, max_iter=1_000_000, max_time=max_time)

base = replace(basename(vtu_path), "_vor.vtu" => "")
outdir = MeshTools.run_outdir("smoothed_$(base)")
MeshTools.save_run_info(outdir, [
    "Script: smooth_mesh_vtu.jl",
    "",
    "Input: $(abspath(vtu_path))",
    "nc = $nc",
    "rtol = $rtol",
    "max_time = $max_time min",
])

label = "mesh_periodic_smoothed_nc$(nc)"
row = MeshTools.save_mesh_level(outdir, smoothed, label, "$nc cells (smoothed)")
MeshTools.print_summary_table([row])
MeshTools.save_summary_csv(joinpath(outdir, "smoothed_metrics_summary.csv"), [row])

println("\nDone! Output written to $(abspath(outdir))/")
