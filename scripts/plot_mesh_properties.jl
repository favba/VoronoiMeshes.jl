# plot_mesh_properties.jl
#
# Reads one or more previously-saved Voronoi meshes (VTU format) and, for each
# one, computes per-cell metrics (see mesh_metrics.jl: area, distortion,
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
# Output: for each mesh "<base>" and each metric in MeshMetrics.METRICS,
# writes "<base>_<metric>.png". Once all meshes are processed, a summary
# table (grid name, resolution and mean/min/max of each metric) is printed
# and saved as "mesh_properties_summary.csv" next to the first processed mesh.
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

include("mesh_metrics.jl")
using .MeshMetrics

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

function save_property_png(filename, mesh, values, title)
    polygons = VoronoiMeshes.create_cell_polygons(mesh)
    fig = Figure(size = (750, 700))
    ax  = Axis(fig[1, 1], title = title, aspect = DataAspect())
    plt = poly!(ax, polygons, color = values, colormap = :viridis, strokewidth = 0.5, strokecolor = (:black, 0.3))
    Colorbar(fig[1, 2], plt)
    hidedecorations!(ax)
    GLMakie.save(filename, fig)
    return nothing
end

# Returns a Dict{Symbol,Any} summary row: :name, :nc, and <metric>_mean/min/max
# for every metric registered in MeshMetrics.METRICS.
function process_mesh(path)
    println("Reading: $path")
    mesh = load_mesh(path)
    nc = length(mesh.cells.position)
    println("  nc = $nc")

    label = base_label(path)
    row = Dict{Symbol, Any}(:name => basename(label), :nc => nc)
    for (mname, mfunc) in MeshMetrics.METRICS
        values = mfunc(mesh)
        m, mn, mx = mean(values), minimum(values), maximum(values)
        println("  $mname: mean = $(round(m, digits=5)), min = $(round(mn, digits=5)), max = $(round(mx, digits=5))")
        row[Symbol(mname, "_mean")] = m
        row[Symbol(mname, "_min")] = mn
        row[Symbol(mname, "_max")] = mx

        filename = "$(label)_$(mname).png"
        save_property_png(filename, mesh, values, "$(basename(label)) — $mname")
        println("  Saved: $filename")
    end
    return row
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

input_paths = if length(ARGS) >= 1
    [ARGS[1]]
else
    manifests = default_manifests()
    isempty(manifests) && error("Usage: julia --project=. plot_mesh_properties.jl <manifest.txt | mesh_vor.vtu>")
    println("No argument given, defaulting to manifests: $(join(manifests, ", "))")
    manifests
end

# Column order: name, nc, then <metric>_mean/min/max for each registered metric.
const SUMMARY_COLUMNS = (
    :name, :nc,
    (Symbol(mname, suffix) for (mname, _) in MeshMetrics.METRICS for suffix in ("_mean", "_min", "_max"))...,
)

fmt_value(val::AbstractString) = val
fmt_value(val::Integer) = string(val)
fmt_value(val::AbstractFloat) = @sprintf("%.5g", val)

function print_summary_table(rows)
    headers = string.(SUMMARY_COLUMNS)
    strs = [[fmt_value(row[key]) for key in SUMMARY_COLUMNS] for row in rows]
    widths = [max(length(headers[i]), maximum(r -> length(r[i]), strs)) for i in eachindex(SUMMARY_COLUMNS)]

    println()
    println(join([rpad(headers[i], widths[i]) for i in eachindex(SUMMARY_COLUMNS)], "  "))
    println(join(["-"^widths[i] for i in eachindex(SUMMARY_COLUMNS)], "  "))
    for s in strs
        println(join([rpad(s[i], widths[i]) for i in eachindex(SUMMARY_COLUMNS)], "  "))
    end
    return nothing
end

function save_summary_csv(filename, rows)
    open(filename, "w") do io
        println(io, join(string.(SUMMARY_COLUMNS), ","))
        for row in rows
            println(io, join((string(row[key]) for key in SUMMARY_COLUMNS), ","))
        end
    end
    return nothing
end

rows = [process_mesh(path) for input_path in input_paths for path in resolve_mesh_paths(input_path)]

if !isempty(rows)
    print_summary_table(rows)
    outdir = dirname(resolve_mesh_paths(input_paths[1])[1])
    csv_path = joinpath(isempty(outdir) ? "." : outdir, "mesh_properties_summary.csv")
    save_summary_csv(csv_path, rows)
    println("\nSaved summary table: $csv_path")
end

println("\nDone!")
