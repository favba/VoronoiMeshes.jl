# Shared build/report pipeline for the build_set_*.jl scripts and
# plot_mesh_properties.jl. Needs a Makie backend (GLMakie/CairoMakie) loaded
# by the including script before `include("mesh_tools.jl")`.

module MeshTools

using VoronoiMeshes: create_cell_polygons, plotmesh!, plotdualmesh!, create_planar_hex_mesh,
                      VoronoiMesh, save_voronoi_to_vtu, save_triangulation_to_vtu, find_obtuse_triangles
using TensorsLiteGeometry: closest
using Statistics: mean
using Dates: now
using Makie: Figure, Axis, Colorbar, poly!, DataAspect, hidedecorations!
using Printf: @sprintf
import Makie

export METRICS, cell_area, cell_area_normalized, cell_distortion, cell_distortion_rms,
       cell_diameter, cell_diameter_normalized, cell_alignment, compute_metrics,
       print_metrics_summary, save_mesh_png, save_property_png, save_all_metric_pngs,
       SUMMARY_COLUMNS, metrics_summary_row, print_summary_table, save_summary_csv,
       save_manifest, rebuild_manifest, numeric_sort_key,
       build_hex_reference, save_mesh_level, finalize_mesh_set,
       center_distance, region_masks, report_region_summary, handle_help, run_outdir,
       save_run_info, obtuse_triangle_count

function handle_help(args, usage)
    if "--help" in args || "-h" in args
        println(usage)
        exit(0)
    end
    return nothing
end

function run_outdir(label)
    outdir = joinpath("output", label)
    mkpath(outdir)
    return outdir
end

function save_run_info(outdir, lines)
    open(joinpath(outdir, "run_info.txt"), "w") do io
        println(io, "Generated: ", now())
        println(io)
        for line in lines
            println(io, line)
        end
    end
    return nothing
end

# Periodic distance from `p` to `center`, correct regardless of where
# `center` sits in the domain (unlike periodic_to_base_point, which only
# wraps `p` itself).
center_distance(p, center, xp, yp) = (v = closest(center, p, xp, yp) - center; hypot(v.x, v.y))

function region_masks(mesh, center, inner_radius, buffer_radius)
    xp, yp = mesh.x_period, mesh.y_period
    d = [center_distance(p, center, xp, yp) for p in mesh.cells.position]
    inner = d .<= inner_radius
    outer = d .> buffer_radius
    buffer = .!inner .& .!outer
    return inner, buffer, outer
end

# dc(x) ∝ ρ(x)^(-1/4) in 2D (Du, Faber & Gunzburger 1999) is only asymptotic,
# so the actual inner/outer dc ratio is measured here rather than assumed.
function report_region_summary(mesh, values, center, inner_radius, buffer_radius, density_ratio)
    inner, buffer, outer = region_masks(mesh, center, inner_radius, buffer_radius)
    dc = values["diameter_dim"]
    n_inner, n_buffer, n_outer = count(inner), count(buffer), count(outer)
    dc_inner = n_inner > 0 ? mean(dc[inner]) : NaN
    dc_outer = n_outer > 0 ? mean(dc[outer]) : NaN
    ratio = dc_inner / dc_outer
    target = density_ratio^(-0.25)
    println("    regions: n_inner=$n_inner, n_buffer=$n_buffer, n_outer=$n_outer")
    println("    dc_inner=$(round(dc_inner, digits=4)), dc_outer=$(round(dc_outer, digits=4)), " *
            "ratio=$(round(ratio, digits=3)) (target $(round(target, digits=3)))")
    if n_buffer > 0
        println("    buffer ring: distortion_rms mean=$(round(mean(values["distortion_rms"][buffer]), digits=4)), " *
                "alignment mean=$(round(mean(values["alignment"][buffer]), digits=4))")
    end
    return (; n_inner, n_buffer, n_outer, dc_inner, dc_outer, ratio)
end

cell_area(mesh) = mesh.cells.area

# Normalized by the mesh's own mean so it's comparable across resolutions
# (raw area/diameter shrink as nc grows).
cell_area_normalized(mesh) = (a = cell_area(mesh); a ./ mean(a))

# (max_edge - min_edge) / mean_edge per cell. Zero for a regular cell.
function cell_distortion(mesh)
    edge_lengths    = mesh.edges.length
    cells_edges     = mesh.cells.edges
    nEdges_per_cell = mesh.cells.nEdges
    nc = length(nEdges_per_cell)

    d = Vector{Float64}(undef, nc)
    for c in 1:nc
        ne   = Int(nEdges_per_cell[c])
        lmin = Inf; lmax = -Inf; lsum = 0.0
        for j in 1:ne
            l    = edge_lengths[cells_edges[c][j]]
            lmin = min(lmin, l)
            lmax = max(lmax, l)
            lsum += l
        end
        d[c] = (lmax - lmin) / (lsum / ne)
    end
    return d
end

# RMS deviation of edge lengths (Peixoto & Barros thesis eq. 2.5, after
# Tomita et al. 2001): S = sqrt(mean((li-lbar)^2))/lbar, lbar = quadratic
# mean of edge lengths. Zero for a regular cell.
function cell_distortion_rms(mesh)
    edge_lengths    = mesh.edges.length
    cells_edges     = mesh.cells.edges
    nEdges_per_cell = mesh.cells.nEdges
    nc = length(nEdges_per_cell)

    S = Vector{Float64}(undef, nc)
    for c in 1:nc
        ne = Int(nEdges_per_cell[c])
        sumsq = 0.0
        for j in 1:ne
            l = edge_lengths[cells_edges[c][j]]
            sumsq += l^2
        end
        lbar = sqrt(sumsq / ne)
        devsq = 0.0
        for j in 1:ne
            l = edge_lengths[cells_edges[c][j]]
            devsq += (l - lbar)^2
        end
        S[c] = sqrt(devsq / ne) / lbar
    end
    return S
end

# Max distance between any two of a cell's vertices, using their periodic
# images closest to the cell center.
function cell_diameter(mesh)
    vert_pos        = mesh.vertices.position
    cell_pos        = mesh.cells.position
    verticesOnCell  = mesh.cells.vertices
    nEdges_per_cell = mesh.cells.nEdges
    xp, yp = mesh.x_period, mesh.y_period
    nc = length(cell_pos)

    diam = Vector{Float64}(undef, nc)
    for c in 1:nc
        ne = Int(nEdges_per_cell[c])
        cp = cell_pos[c]
        vs = verticesOnCell[c]
        local_vertices = ntuple(j -> closest(cp, vert_pos[vs[j]], xp, yp), ne)
        dmax = 0.0
        for i in 1:ne, j in (i + 1):ne
            dv = local_vertices[i] - local_vertices[j]
            dmax = max(dmax, sqrt(dv.x^2 + dv.y^2))
        end
        diam[c] = dmax
    end
    return diam
end

cell_diameter_normalized(mesh) = (d = cell_diameter(mesh); d ./ mean(d))

# Alignment index Ξ (Peixoto & Barros 2013, Prop. 3.1.5): 0 for a cell whose
# opposite edges are parallel and equal length, growing as it departs from
# that symmetry. Odd-sided cells (no even-alignment notion) get Ξ = 0.
function cell_alignment(mesh)
    vert_pos        = mesh.vertices.position
    cell_pos        = mesh.cells.position
    verticesOnCell  = mesh.cells.vertices
    nEdges_per_cell = mesh.cells.nEdges
    xp, yp = mesh.x_period, mesh.y_period
    nc = length(cell_pos)

    align = Vector{Float64}(undef, nc)
    for c in 1:nc
        n = Int(nEdges_per_cell[c])
        if isodd(n)
            align[c] = 0.0
            continue
        end
        cp = cell_pos[c]
        vs = verticesOnCell[c]
        P = ntuple(j -> closest(cp, vert_pos[vs[j]], xp, yp), n)

        dist(i, j) = (a = P[mod1(i, n)]; b = P[mod1(j, n)]; hypot(a.x - b.x, a.y - b.y))

        dbar = sum(dist(i, i + 1) for i in 1:n) / n
        half = n ÷ 2
        s = 0.0
        for i in 1:half
            s += abs(dist(i + 1 + half, i) - dist(i + half, i + 1))
            s += abs(dist(i + 1, i) - dist(i + half + 1, i + half))
        end
        align[c] = s / (dbar * n)
    end
    return align
end

const METRICS = (
    ("diameter_dim", cell_diameter),
    ("area", cell_area_normalized),
    ("distortion", cell_distortion),
    ("distortion_rms", cell_distortion_rms),
    ("diameter", cell_diameter_normalized),
    ("alignment", cell_alignment),
)

compute_metrics(mesh) = Dict(mname => mfunc(mesh) for (mname, mfunc) in METRICS)

obtuse_triangle_count(mesh) = length(find_obtuse_triangles(mesh)), length(mesh.vertices.cells)

function print_metrics_summary(mesh, values, label)
    nc = length(mesh.cells.position)
    n_obtuse, n_triangles = obtuse_triangle_count(mesh)
    println("  Metrics ($label): nc = $nc, x_period = $(mesh.x_period), y_period = $(mesh.y_period)")
    println("    obtuse triangles: $n_obtuse / $n_triangles")
    for (mname, _) in METRICS
        v = values[mname]
        println("    $mname: mean = $(round(mean(v), digits=5)), min = $(round(minimum(v), digits=5)), max = $(round(maximum(v), digits=5))")
    end
    return nothing
end

function save_mesh_png(filename, mesh, label)
    fig = Figure(size=(700, 700))
    ax = Axis(fig[1, 1], title=label, aspect=DataAspect())
    plotdualmesh!(ax, mesh)
    plotmesh!(ax, mesh)
    hidedecorations!(ax)
    Makie.save(filename, fig)
    return nothing
end

function save_property_png(filename, mesh, values, title)
    polygons = create_cell_polygons(mesh)
    fig = Figure(size=(750, 700))
    ax = Axis(fig[1, 1], title=title, aspect=DataAspect())
    plt = poly!(ax, polygons, color=values, colormap=:viridis, strokewidth=0.5, strokecolor=(:black, 0.3))
    Colorbar(fig[1, 2], plt)
    hidedecorations!(ax)
    Makie.save(filename, fig)
    return nothing
end

function save_all_metric_pngs(label, mesh, values)
    for (mname, _) in METRICS
        filename = "$(label)_$(mname).png"
        save_property_png(filename, mesh, values[mname], "$(basename(label)) — $mname")
        println("  Saved: $filename")
    end
    return nothing
end

const SUMMARY_COLUMNS = (
    :name, :nc, :n_obtuse, :n_triangles,
    (Symbol(mname, suffix) for (mname, _) in METRICS for suffix in ("_mean", "_min", "_max"))...,
)

function metrics_summary_row(mesh, values, name)
    nc = length(mesh.cells.position)
    n_obtuse, n_triangles = obtuse_triangle_count(mesh)
    row = Dict{Symbol, Any}(:name => name, :nc => nc, :n_obtuse => n_obtuse, :n_triangles => n_triangles)
    for (mname, _) in METRICS
        v = values[mname]
        row[Symbol(mname, "_mean")] = mean(v)
        row[Symbol(mname, "_min")] = minimum(v)
        row[Symbol(mname, "_max")] = maximum(v)
    end
    return row
end

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

# Appends rather than overwrites, so a later run adding more levels/scales
# doesn't wipe out rows already on disk.
function save_summary_csv(filename, rows)
    write_header = !isfile(filename)
    open(filename, write_header ? "w" : "a") do io
        write_header && println(io, join(string.(SUMMARY_COLUMNS), ","))
        for row in rows
            println(io, join((string(row[key]) for key in SUMMARY_COLUMNS), ","))
        end
    end
    return nothing
end

function save_manifest(filename, vtu_names)
    open(filename, "w") do io
        for name in vtu_names
            println(io, name)
        end
    end
    return nothing
end

# Scans `outdir` for files matching `pattern` (not just what this run
# produced, so leftovers from earlier runs are picked up too) and appends
# any not already listed to the manifest, ordered by `sort_key`.
function rebuild_manifest(outdir, manifest_name, pattern, sort_key)
    files = filter(f -> occursin(pattern, f), readdir(outdir))
    isempty(files) && return 0
    sort!(files, by=sort_key)

    manifest_path = joinpath(outdir, manifest_name)
    if !isfile(manifest_path)
        save_manifest(manifest_path, files)
        return length(files)
    end

    existing = Set(readlines(manifest_path))
    new_files = filter(f -> !(f in existing), files)
    if !isempty(new_files)
        open(manifest_path, "a") do io
            for f in new_files
                println(io, f)
            end
        end
    end
    return length(existing) + length(new_files)
end

numeric_sort_key(pattern) = f -> Tuple(parse(Float64, g) for g in match(pattern, f).captures)

# create_planar_hex_mesh rounds nc to fit an integer grid, so its own domain
# isn't exactly xperiod x yperiod; reusing its generator count against the
# exact domain lets Lloyd relaxation fill it precisely.
function build_hex_reference(nc, xperiod, yperiod; density=nothing, kwargs...)
    dc = sqrt(xperiod * yperiod / nc)
    hex_mesh = create_planar_hex_mesh(xperiod, yperiod, dc)
    mesh = density === nothing ?
        VoronoiMesh(hex_mesh.cells.position, xperiod, yperiod; kwargs...) :
        VoronoiMesh(hex_mesh.cells.position, xperiod, yperiod; density, kwargs...)
    return mesh, dc
end

# `region`, if given as (center, inner_radius, buffer_radius, density_ratio),
# also reports per-region stats on the same precomputed metrics.
function save_mesh_level(outdir, mesh, label, title=label; region=nothing)
    save_voronoi_to_vtu(joinpath(outdir, "$(label)_vor.vtu"), mesh)
    save_triangulation_to_vtu(joinpath(outdir, "$(label)_tri.vtu"), mesh)
    save_mesh_png(joinpath(outdir, "$(label).png"), mesh, title)
    println("  Saved: $label")
    values = compute_metrics(mesh)
    print_metrics_summary(mesh, values, title)
    region === nothing || report_region_summary(mesh, values, region...)
    save_all_metric_pngs(joinpath(outdir, label), mesh, values)
    return metrics_summary_row(mesh, values, label)
end

function finalize_mesh_set(outdir, kind, rows, pattern, sort_key)
    if isempty(rows)
        println("\nNo meshes were built this run; skipping summary table/CSV.")
    else
        print_summary_table(rows)
        csv_path = joinpath(outdir, "$(kind)_metrics_summary.csv")
        save_summary_csv(csv_path, rows)
        println("\nSaved summary table: $csv_path")
    end

    manifest_name = "$(kind)_voronoi_meshes.txt"
    manifest_path = joinpath(outdir, manifest_name)
    n = rebuild_manifest(outdir, manifest_name, pattern, sort_key)
    if n > 0
        println("Saved manifest: $manifest_path ($n meshes)")
        n != length(rows) && println(
            "  Note: manifest covers $n mesh(es) found on disk in $outdir; " *
            "this run's summary/CSV covers $(length(rows)) of them.",
        )
    else
        println("No $kind meshes found in $outdir matching the expected naming pattern; " *
                "manifest at $manifest_path left unchanged.")
    end

    println("\nDone! Output written to $(abspath(outdir))/")
    return nothing
end

end # module
