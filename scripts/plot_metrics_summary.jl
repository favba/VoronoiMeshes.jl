# Reads a "<kind>_metrics_summary.csv" (as written by mesh_tools.jl's
# save_summary_csv) and plots distortion_rms/alignment/area/diameter against
# either cell count (nc varies row to row) or perturbation strength d (nc
# fixed, "name" column encodes "..._d<strength>") — auto-detected.
#
# Usage:
#   julia --project=. plot_metrics_summary.jl <metrics_summary.csv>
#
# Output: "<csv_basename>_convergence.pdf"/".eps" next to the input CSV.

using CairoMakie
using Printf

const NUMERIC_COLUMNS = (:nc, :diameter_dim_mean, :distortion_rms_mean, :alignment_mean, :area_min, :area_max, :diameter_min, :diameter_max)

function read_summary_csv(path)
    lines = readlines(path)
    header = split(lines[1], ",")
    col_idx = Dict(Symbol(h) => i for (i, h) in enumerate(header))
    haskey(col_idx, :name) || error("Column 'name' not found in $path")
    for c in NUMERIC_COLUMNS
        haskey(col_idx, c) || error("Column '$c' not found in $path")
    end

    rows = [split(line, ",") for line in lines[2:end] if !isempty(strip(line))]
    data = Dict{Symbol, Any}(c => [parse(Float64, r[col_idx[c]]) for r in rows] for c in NUMERIC_COLUMNS)
    data[:name] = [String(r[col_idx[:name]]) for r in rows]

    order = sortperm(data[:nc])
    for c in (NUMERIC_COLUMNS..., :name)
        data[c] = data[c][order]
    end
    return data
end

# Extracts the perturbation/distortion strength from a mesh name of the form
# "..._d<strength>" (the naming convention build_set_global_distortion_meshes.jl
# uses for each level), e.g. "mesh_periodic_global_distortion_nc80_d0.45" -> 0.45.
function perturbation_strength(name)
    m = match(r"_d([0-9.]+)$", name)
    m === nothing && error("Could not parse perturbation strength 'd' from name: $name")
    return parse(Float64, m.captures[1])
end

# Okabe-Ito colorblind-safe categorical colors, one per plotted series.
const COLOR_DISTORTION = :black
const COLOR_ALIGNMENT  = RGBf(0.90, 0.62, 0.0)
const COLOR_AREA       = RGBf(0.0, 0.45, 0.70)
const COLOR_DIAMETER   = RGBf(0.0, 0.62, 0.45)

# Domain side length 1 is treated as a great-circle circumference of Earth.
const EARTH_RADIUS_KM = 6371.0
const DOMAIN_UNIT_TO_KM = 2 * pi * EARTH_RADIUS_KM

function plot_metric_series!(fig, x, data)
    ax = Axis(fig[2, 1], ylabel = "Normalized metric value")

    area_ratio = data[:area_min] ./ data[:area_max]
    diameter_ratio = data[:diameter_min] ./ data[:diameter_max]

    lines!(ax, x, data[:distortion_rms_mean], color = COLOR_DISTORTION, linewidth = 2)
    scatter!(ax, x, data[:distortion_rms_mean], color = COLOR_DISTORTION, markersize = 9, label = "distortion")

    lines!(ax, x, data[:alignment_mean], color = COLOR_ALIGNMENT, linewidth = 2)
    scatter!(ax, x, data[:alignment_mean], color = COLOR_ALIGNMENT, markersize = 9, label = "misalignment")

    lines!(ax, x, area_ratio, color = COLOR_AREA, linewidth = 2)
    scatter!(ax, x, area_ratio, color = COLOR_AREA, markersize = 9, label = "area (min/max)")

    lines!(ax, x, diameter_ratio, color = COLOR_DIAMETER, linewidth = 2)
    scatter!(ax, x, diameter_ratio, color = COLOR_DIAMETER, markersize = 9, label = "diameter (min/max)")

    axislegend(ax, position = :lt, framevisible = true, backgroundcolor = (:white, 0.75))
    return ax
end

function plot_resolution_sweep(fig, data)
    nc = data[:nc]
    ax = plot_metric_series!(fig, nc, data)
    ax.xlabel = "Number of cells (resolution)"
    ax.xscale = log2
    ax.xticks = (nc, string.(round.(Int, nc)))

    ax_top = Axis(fig[2, 1],
        xaxisposition = :top, xscale = log2,
        xlabel = "Mean cell diameter (km)",
        xticks = (nc, [@sprintf("%.0f", d * DOMAIN_UNIT_TO_KM) for d in data[:diameter_dim_mean]]))
    hideydecorations!(ax_top)
    hidespines!(ax_top, :l, :r, :b)
    ax_top.xgridvisible = false
    linkxaxes!(ax, ax_top)
    return nothing
end

function plot_perturbation_sweep(fig, data)
    nc = round(Int, data[:nc][1])
    d = perturbation_strength.(data[:name])

    order = sortperm(d)
    d = d[order]
    sorted_data = Dict(c => data[c][order] for c in keys(data))

    ax = plot_metric_series!(fig, d, sorted_data)
    ax.xlabel = "Perturbation strength d (fixed resolution of $nc cells)"
    return nothing
end

function plot_convergence(csv_path)
    data = read_summary_csv(csv_path)

    fig = Figure(size = (800, 620))
    Label(fig[1, 1], basename(csv_path), fontsize = 16, font = :bold, tellwidth = false)

    if length(unique(data[:nc])) == 1
        plot_perturbation_sweep(fig, data)
    else
        plot_resolution_sweep(fig, data)
    end

    base_path = replace(csv_path, r"\.csv$" => "_convergence")
    for ext in ("pdf", "eps")
        out_path = "$(base_path).$(ext)"
        save(out_path, fig)
        println("Saved: $out_path")
    end
    return base_path
end

length(ARGS) == 1 || error("Usage: julia --project=. plot_metrics_summary.jl <metrics_summary.csv>")
plot_convergence(ARGS[1])
