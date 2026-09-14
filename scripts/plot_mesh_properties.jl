# Computes per-cell metrics for previously-saved meshes and plots each as a
# colored Voronoi diagram PNG.
#
# Usage:
#   julia --project=. plot_mesh_properties.jl [manifest.txt | mesh_vor.vtu]
#
# With no argument, processes every "*_voronoi_meshes.txt" manifest found
# under output/. A manifest lists one "*_vor.vtu" path per line, relative to
# its own directory; a bare mesh path works too ("<base>_vor.vtu" or
# "<base>.vtu"). Writes "<base>_<metric>.png" per mesh, plus a summary CSV
# next to each input ("<kind>_metrics_summary.csv" for a manifest,
# "mesh_properties_summary.csv" for a bare mesh).

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

function default_manifests()
    isdir("output") || return String[]
    found = String[]
    for (dir, _, files) in walkdir("output")
        for f in files
            endswith(f, "_voronoi_meshes.txt") && push!(found, joinpath(dir, f))
        end
    end
    sort!(found)
    return found
end

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
