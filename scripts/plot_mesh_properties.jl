# plot_mesh_properties.jl
#
# Reads one or more previously-saved Voronoi meshes (VTU format) and, for each
# one, computes per-cell metrics (see mesh_tools.jl: area, distortion,
# diameter, ...), then plots each property as a colored Voronoi diagram saved
# to a PNG.
#
# Usage:
#   julia --project=. plot_mesh_properties.jl [manifest.txt | mesh_vor.vtu]
#
#   With no argument, processes every "output/*_voronoi_meshes.txt" manifest
#   produced by the build_set_*.jl scripts.
#
#   - manifest.txt : a text file listing one "*_vor.vtu" path per line
#                    (blank lines and lines starting with '#' are ignored;
#                    relative paths are resolved against the manifest's
#                    own directory). This is the format written by the
#                    build_set_*.jl scripts into scripts/output/*.txt.
#   - mesh_vor.vtu : a single mesh, given either as "<base>_vor.vtu" (the
#                    matching "<base>_tri.vtu" is derived automatically) or
#                    as a bare "<base>.vtu" base name.
#
# Output: for each mesh "<base>" and each metric in MeshTools.METRICS,
# writes "<base>_<metric>.png". Each input (one manifest, or one bare mesh) is
# summarized separately: a table is printed and a CSV is saved next to the
# processed mesh(es), named after the manifest — "<kind>_voronoi_meshes.txt"
# produces "<kind>_metrics_summary.csv" (matching the name the build_set_*.jl
# scripts themselves use); a bare mesh path produces "mesh_properties_summary.csv".
#
# Note: PNG export uses GLMakie. On headless servers swap to CairoMakie
# (add it to the environment and replace `using GLMakie` below).

using VoronoiMeshes
using TensorsLite
using TensorsLiteGeometry
using ReadVTK
using GLMakie
using Statistics
using Printf

include("mesh_tools.jl")
using .MeshTools

function load_mesh(path)
    if endswith(path, "_vor.vtu")
        tri_path = replace(path, "_vor.vtu" => "_tri.vtu")
        return VoronoiMeshes.read_from_vtu(path, tri_path)
    else
        return VoronoiMeshes.read_from_vtu(path)
    end
end

function base_label(path)
    name = basename(path)
    name = replace(name, "_vor.vtu" => "", ".vtu" => "")
    return joinpath(dirname(path), name)
end

# Returns a Dict{Symbol,Any} summary row: :name, :nc, and <metric>_mean/min/max
# for every metric registered in MeshTools.METRICS.
function process_mesh(path)
    println("Reading: $path")
    mesh = load_mesh(path)
    label = base_label(path)

    values = MeshTools.compute_metrics(mesh)
    MeshTools.print_metrics_summary(mesh, values, basename(label))
    MeshTools.save_all_metric_pngs(label, mesh, values)

    return MeshTools.metrics_summary_row(mesh, values, basename(label))
end

function resolve_mesh_paths(input_path)
    if endswith(input_path, ".txt")
        dir = dirname(input_path)
        lines = readlines(input_path)
        paths = String[]
        for line in lines
            l = strip(line)
            (isempty(l) || startswith(l, "#")) && continue
            push!(paths, isabspath(l) ? l : joinpath(dir, l))
        end
        return paths
    else
        return [input_path]
    end
end

# Default: no argument given -> process every manifest already in output/.
default_manifests() = filter(isfile, joinpath.("output", (
    "regular_voronoi_meshes.txt", "irregular_voronoi_meshes.txt", "refined_voronoi_meshes.txt",
)))

# Derives this input's own summary CSV name: "<kind>_voronoi_meshes.txt" (the
# convention the build_set_*.jl scripts write) -> "<kind>_metrics_summary.csv";
# any other ".txt" -> "<name>_summary.csv"; a bare mesh path -> the generic name.
function summary_csv_name(input_path)
    endswith(input_path, ".txt") || return "mesh_properties_summary.csv"
    base = replace(basename(input_path), "_voronoi_meshes.txt" => "", ".txt" => "")
    return "$(base)_metrics_summary.csv"
end

input_paths = if length(ARGS) >= 1
    [ARGS[1]]
else
    manifests = default_manifests()
    isempty(manifests) && error("Usage: julia --project=. plot_mesh_properties.jl <manifest.txt | mesh_vor.vtu>")
    println("No argument given, defaulting to manifests: $(join(manifests, ", "))")
    manifests
end

for input_path in input_paths
    paths = resolve_mesh_paths(input_path)
    rows = [process_mesh(path) for path in paths]
    isempty(rows) && continue

    MeshTools.print_summary_table(rows)
    outdir = dirname(paths[1])
    csv_path = joinpath(isempty(outdir) ? "." : outdir, summary_csv_name(input_path))
    MeshTools.save_summary_csv(csv_path, rows)
    println("\nSaved summary table: $csv_path")
end

println("\nDone!")
