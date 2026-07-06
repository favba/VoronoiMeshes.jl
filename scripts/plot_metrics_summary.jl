# plot_metrics_summary.jl
#
# Reads one "<kind>_metrics_summary.csv" file (as produced by the
# build_set_*.jl scripts / mesh_tools.jl's save_summary_csv, one row per
# mesh resolution) and plots a paper-ready convergence figure: number of
# cells (resolution) on the bottom x-axis (with a linked top x-axis relabeled
# in the corresponding dimensional mean cell diameter), and a set of
# normalized quality metrics on a shared y-axis — distortion_rms, alignment,
# and the min/max ratio (min/max, a measure of spread) of area and diameter.
#
# Usage:
#   julia --project=. plot_metrics_summary.jl <metrics_summary.csv>
#
# Output: "<csv_basename>_convergence.pdf" next to the input CSV.
#
# Uses CairoMakie (vector output) rather than GLMakie, since this produces a
# static publication figure rather than an interactive/rendered mesh view.

using CairoMakie
using Printf

const COLUMNS = (:nc, :diameter_dim_mean, :distortion_rms_mean, :alignment_mean, :area_min, :area_max, :diameter_min, :diameter_max)

# Reads `path` (the plain, unquoted CSV format written by
# MeshTools.save_summary_csv) into a Dict{Symbol,Vector{Float64}} for the
# columns in COLUMNS, sorted by ascending nc.
function read_summary_csv(path)
    lines = readlines(path)
    header = split(lines[1], ",")
    col_idx = Dict(Symbol(h) => i for (i, h) in enumerate(header))
    for c in COLUMNS
        haskey(col_idx, c) || error("Column '$c' not found in $path")
    end

    rows = [split(line, ",") for line in lines[2:end] if !isempty(strip(line))]
    data = Dict(c => [parse(Float64, r[col_idx[c]]) for r in rows] for c in COLUMNS)

    order = sortperm(data[:nc])
    for c in COLUMNS
        data[c] = data[c][order]
    end
    return data
end

# Okabe-Ito colorblind-safe categorical colors, one per plotted series.
const COLOR_DISTORTION = :black
const COLOR_ALIGNMENT  = RGBf(0.90, 0.62, 0.0)
const COLOR_AREA       = RGBf(0.0, 0.45, 0.70)
const COLOR_DIAMETER   = RGBf(0.0, 0.62, 0.45)

# Mean earth radius (km); the periodic planar domain (side length 1) is
# treated as a great-circle circumference of the sphere these meshes will
# eventually be mapped onto, so a domain-unit length converts to km via
# 2*pi*EARTH_RADIUS_KM.
const EARTH_RADIUS_KM = 6371.0
const DOMAIN_UNIT_TO_KM = 2 * pi * EARTH_RADIUS_KM

function plot_convergence(csv_path)
    data = read_summary_csv(csv_path)
    nc = data[:nc]
    area_ratio = data[:area_min] ./ data[:area_max]
    diameter_ratio = data[:diameter_min] ./ data[:diameter_max]

    fig = Figure(size = (800, 620))
    Label(fig[1, 1], basename(csv_path), fontsize = 16, font = :bold, tellwidth = false)

    ax = Axis(fig[2, 1],
        xlabel = "Number of cells (resolution)", ylabel = "Normalized metric value",
        xscale = log2, xticks = (nc, string.(round.(Int, nc))))

    # Top twin x-axis: same tick positions as nc, relabeled with the
    # corresponding dimensional cell diameter — not a second series, just the
    # existing x-axis expressed in physical units alongside the cell count.
    ax_top = Axis(fig[2, 1],
        xaxisposition = :top, xscale = log2,
        xlabel = "Mean cell diameter (km)",
        xticks = (nc, [@sprintf("%.3g", d * DOMAIN_UNIT_TO_KM) for d in data[:diameter_dim_mean]]))
    hideydecorations!(ax_top)
    hidespines!(ax_top, :l, :r, :b)
    ax_top.xgridvisible = false
    linkxaxes!(ax, ax_top)

    lines!(ax, nc, data[:distortion_rms_mean], color = COLOR_DISTORTION, linewidth = 2)
    scatter!(ax, nc, data[:distortion_rms_mean], color = COLOR_DISTORTION, markersize = 9, label = "distortion")

    lines!(ax, nc, data[:alignment_mean], color = COLOR_ALIGNMENT, linewidth = 2)
    scatter!(ax, nc, data[:alignment_mean], color = COLOR_ALIGNMENT, markersize = 9, label = "misalignment")

    lines!(ax, nc, area_ratio, color = COLOR_AREA, linewidth = 2)
    scatter!(ax, nc, area_ratio, color = COLOR_AREA, markersize = 9, label = "area (min/max)")

    lines!(ax, nc, diameter_ratio, color = COLOR_DIAMETER, linewidth = 2)
    scatter!(ax, nc, diameter_ratio, color = COLOR_DIAMETER, markersize = 9, label = "diameter (min/max)")

    axislegend(ax, position = :lt, framevisible = true, backgroundcolor = (:white, 0.75))

    out_path = replace(csv_path, r"\.csv$" => "_convergence.pdf")
    save(out_path, fig)
    println("Saved: $out_path")
    return out_path
end

length(ARGS) == 1 || error("Usage: julia --project=. plot_metrics_summary.jl <metrics_summary.csv>")
plot_convergence(ARGS[1])
