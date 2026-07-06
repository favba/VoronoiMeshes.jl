# mesh_tools.jl
#
# Shared build/report pipeline for the build_set_*.jl scripts and
# plot_mesh_properties.jl: per-cell quality metrics, the plotting helpers that
# render them, and the printing/CSV/manifest bookkeeping every mesh-building
# script needs.
#
# Each metric function takes a VoronoiMesh and returns a Vector with one
# value per cell. Add new metrics here and register them in `METRICS` so
# they automatically show up in scripts that iterate over it.
#
# Note: the plotting helpers below need a Makie backend (GLMakie/CairoMakie)
# already loaded by the including script before `include("mesh_tools.jl")`.

module MeshTools

using VoronoiMeshes: create_cell_polygons, plotmesh!, plotdualmesh!, create_planar_hex_mesh,
                      VoronoiMesh, save_voronoi_to_vtu, save_triangulation_to_vtu
using TensorsLiteGeometry: closest
using Statistics: mean
using Makie: Figure, Axis, Colorbar, poly!, DataAspect, hidedecorations!
using Printf: @sprintf
import Makie

export METRICS, cell_area, cell_area_normalized, cell_distortion, cell_distortion_rms,
       cell_diameter, cell_diameter_normalized, cell_alignment, compute_metrics,
       print_metrics_summary, save_mesh_png, save_property_png, save_all_metric_pngs,
       SUMMARY_COLUMNS, metrics_summary_row, print_summary_table, save_summary_csv,
       save_manifest, rebuild_manifest, numeric_sort_key,
       build_hex_reference, save_mesh_level, finalize_mesh_set

# Cell area, already computed and cached by VoronoiMeshes.
cell_area(mesh) = mesh.cells.area

# Cell area relative to the mesh's own mean area. Raw area/diameter shrink with
# resolution (area ~ 1/nc, diameter ~ 1/sqrt(nc)), so comparing their raw values
# across different grid resolutions is meaningless. Normalizing by the mesh's
# own mean turns them into a dimensionless "relative size" (1.0 = average cell),
# which is what's actually comparable across resolutions. distortion/distortion_rms/
# alignment are already dimensionless ratios and need no such normalization.
cell_area_normalized(mesh) = (a = cell_area(mesh); a ./ mean(a))

# Per-cell irregularity: (max_edge - min_edge) / mean_edge over each cell's edges.
# Zero for a perfectly regular cell, grows with distortion.
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

# Per-cell distortion S, as formally defined in Peixoto & Barros thesis, eq. (2.5)
# (following Tomita et al. 2001): S = sqrt(mean((li - lbar)^2)) / lbar, where
# lbar = sqrt(mean(li^2)) is the quadratic mean of the cell's edge lengths.
# This is an RMS deviation of edge lengths, distinct from the simpler
# (max-min)/mean ratio used by `cell_distortion` above. Zero for a regular cell.
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

# Per-cell diameter: max distance between any two of the cell's vertices,
# using their periodic images closest to the cell center.
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

# Cell diameter relative to the mesh's own mean diameter (see cell_area_normalized).
cell_diameter_normalized(mesh) = (d = cell_diameter(mesh); d ./ mean(d))

# Per-cell alignment index Ξ (Peixoto & Barros 2013, Prop. 3.1.5 / thesis eq. 3.1),
# adapted from the sphere to this package's periodic planar meshes.
# A polygon with an even number of vertices is "aligned" (opposite edges parallel
# and of equal length) iff Ξ = 0; Ξ grows as the cell departs from that symmetry.
# Cells with an odd number of edges (e.g. pentagons/heptagons from local defects)
# have no even-alignment notion, so Ξ is set to 0 for them, matching the thesis'
# convention of assigning pentagons a null index purely so they can be plotted
# alongside hexagons.
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

# Ordered (name, function) pairs used by scripts to compute and report
# every available per-cell metric without hardcoding the list.
const METRICS = (
    ("area", cell_area_normalized),
    ("distortion", cell_distortion),
    ("distortion_rms", cell_distortion_rms),
    ("diameter", cell_diameter_normalized),
    ("alignment", cell_alignment),
)

# Computes every registered metric for `mesh` once, returning a
# Dict{String,Vector} keyed by METRICS name. Callers that need more than one
# of print_metrics_summary/save_all_metric_pngs/metrics_summary_row should
# compute this once and pass it to all of them, rather than letting each one
# recompute every metric (some, like cell_diameter/cell_alignment, are not
# O(1) per cell) from scratch.
compute_metrics(mesh) = Dict(mname => mfunc(mesh) for (mname, mfunc) in METRICS)

# Prints mean/min/max of every registered metric, prefixed by `label`, from
# precomputed `values` (see `compute_metrics`).
function print_metrics_summary(mesh, values, label)
    nc = length(mesh.cells.position)
    println("  Metrics ($label): nc = $nc, x_period = $(mesh.x_period), y_period = $(mesh.y_period)")
    for (mname, _) in METRICS
        v = values[mname]
        println("    $mname: mean = $(round(mean(v), digits=5)), min = $(round(minimum(v), digits=5)), max = $(round(maximum(v), digits=5))")
    end
    return nothing
end

# Plain overlay plot: Voronoi diagram (blue) over its dual triangulation (orange).
function save_mesh_png(filename, mesh, label)
    fig = Figure(size=(700, 700))
    ax = Axis(fig[1, 1], title=label, aspect=DataAspect())
    plotdualmesh!(ax, mesh)
    plotmesh!(ax, mesh)
    hidedecorations!(ax)
    Makie.save(filename, fig)
    return nothing
end

# Voronoi diagram colored by a per-cell scalar (e.g. one of the METRICS).
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

# Saves "<label>_<metric>.png" for every registered metric, from precomputed
# `values` (see `compute_metrics`).
function save_all_metric_pngs(label, mesh, values)
    for (mname, _) in METRICS
        filename = "$(label)_$(mname).png"
        save_property_png(filename, mesh, values[mname], "$(basename(label)) — $mname")
        println("  Saved: $filename")
    end
    return nothing
end

# Column order shared by print_summary_table/save_summary_csv: name, nc, then
# <metric>_mean/min/max for each metric registered in METRICS.
const SUMMARY_COLUMNS = (
    :name, :nc,
    (Symbol(mname, suffix) for (mname, _) in METRICS for suffix in ("_mean", "_min", "_max"))...,
)

# One summary row (Dict{Symbol,Any}) for `mesh`, identified by `name` in the
# resulting table/CSV (e.g. a mesh label or level tag), from precomputed
# `values` (see `compute_metrics`).
function metrics_summary_row(mesh, values, name)
    nc = length(mesh.cells.position)
    row = Dict{Symbol, Any}(:name => name, :nc => nc)
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

# `filename` is caller-supplied (not hardcoded here) so each script can pick a
# name that won't collide with the other build_set_*.jl / plot_mesh_properties.jl
# summaries sharing the same output directory.
function save_summary_csv(filename, rows)
    open(filename, "w") do io
        println(io, join(string.(SUMMARY_COLUMNS), ","))
        for row in rows
            println(io, join((string(row[key]) for key in SUMMARY_COLUMNS), ","))
        end
    end
    return nothing
end

# Writes one "<mesh>_vor.vtu" name per line — a manifest that plot_mesh_properties.jl
# reads (as "output/*_voronoi_meshes.txt") to process every mesh a build_set_*.jl
# script produced, without having to guess at wildcards. `vtu_names` are written
# as given, resolved relative to the manifest's own directory by the reader.
function save_manifest(filename, vtu_names)
    open(filename, "w") do io
        for name in vtu_names
            println(io, name)
        end
    end
    return nothing
end

# Scans `outdir` for "*_vor.vtu" filenames matching `pattern` and (re)writes
# "<outdir>/<manifest_name>" listing all of them, ordered by `sort_key`
# (typically the level/scale index parsed out of the filename, e.g.
# `f -> parse(Int, match(pattern, f)[1])`).
#
# Rebuilding from whatever mesh files are actually on disk — rather than only
# the ones a single run just produced — means the manifest stays complete even
# when it should include meshes left over from an earlier run (different nc,
# different random seed, or built before manifest-writing existed).
function rebuild_manifest(outdir, manifest_name, pattern, sort_key)
    files = filter(f -> occursin(pattern, f), readdir(outdir))
    isempty(files) && return 0
    sort!(files, by=sort_key)
    save_manifest(joinpath(outdir, manifest_name), files)
    return length(files)
end

# Builds a sort_key function for rebuild_manifest/finalize_mesh_set: matches
# `pattern` against a filename and returns its captured groups, parsed as
# Ints, as a tuple (in capture-group order). Centralizes the "regex match +
# parse capture groups" idiom so each build script only needs to supply its
# own naming pattern, not re-derive a matching closure by hand.
numeric_sort_key(pattern) = f -> Tuple(parse(Int, g) for g in match(pattern, f).captures)

# Builds the regular hex-mesh reference shared by the regular and irregular
# build scripts (their common "Level 0"): a hex mesh at ~nc cells, rebuilt
# against the exact periodic domain. create_planar_hex_mesh rounds the cell
# count to fit an integer number of hex rows/columns, so its own returned
# domain isn't exactly xperiod x yperiod; reusing its generator count (same
# cell count) against the exact domain lets Lloyd relaxation (VoronoiMesh's
# default) spread the generators to fill it precisely.
function build_hex_reference(nc, xperiod, yperiod)
    dc = sqrt(xperiod * yperiod / nc)
    hex_mesh = create_planar_hex_mesh(xperiod, yperiod, dc)
    mesh = VoronoiMesh(hex_mesh.cells.position, xperiod, yperiod)
    return mesh, dc
end

# Saves one mesh "level"/scale to disk — Voronoi + triangulation VTU, a plain
# overlay PNG, per-metric colored PNGs — and prints the metrics summary. This
# is the common bundle of outputs every build_set_*.jl script produces per
# mesh it generates. Returns the metrics_summary_row for `label`, which
# callers accumulate into `rows` for the final summary table/CSV (see
# `finalize_mesh_set`).
function save_mesh_level(outdir, mesh, label, title=label)
    save_voronoi_to_vtu(joinpath(outdir, "$(label)_vor.vtu"), mesh)
    save_triangulation_to_vtu(joinpath(outdir, "$(label)_tri.vtu"), mesh)
    save_mesh_png(joinpath(outdir, "$(label).png"), mesh, title)
    println("  Saved: $label")
    values = compute_metrics(mesh)
    print_metrics_summary(mesh, values, title)
    save_all_metric_pngs(joinpath(outdir, label), mesh, values)
    return metrics_summary_row(mesh, values, label)
end

# Wraps up a mesh_set build: prints/saves the combined summary table + CSV for
# `rows`, and (re)writes the `kind`'s manifest by rescanning `outdir` for every
# matching mesh (see `rebuild_manifest`) rather than just the ones this run
# produced. `kind` names the output files: "<kind>_metrics_summary.csv" and
# "<kind>_voronoi_meshes.txt". The common "wrap up a build" step shared by all
# three build_set_*.jl scripts.
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
