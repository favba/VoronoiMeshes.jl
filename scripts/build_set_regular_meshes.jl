# build_set_regular_meshes.jl
#
# Builds a series of regular Voronoi meshes by successive refinement starting
# from a regular hexagonal mesh.
#
#   Level 0 : regular hex mesh from create_planar_hex_mesh (~nc_ref cells)
#   Level i : uses cell centres + edge midpoints of level i-1 as generators,
#             then runs Lloyd relaxation to convergence.
#             This roughly quadruples the cell count at each level.
#
# Usage:
#   julia --project=. build_set_regular_meshes.jl [nc_ref] [num_levels]
#
# Defaults: nc_ref=16, num_levels=4

using VoronoiMeshes
using DelaunayTriangulation
using TensorsLite
using WriteVTK
using GLMakie
using Statistics

const X_PERIOD = 1.0
const Y_PERIOD = 1.0

# Per-cell irregularity: (max_edge - min_edge) / mean_edge, averaged over all cells.
# Zero for perfectly regular cells, grows with distortion.
function distortion_metric(mesh)
    edge_lengths    = mesh.edges.length
    cells_edges     = mesh.cells.edges
    nEdges_per_cell = mesh.cells.nEdges
    nc = length(nEdges_per_cell)

    cell_d = Vector{Float64}(undef, nc)
    for c in 1:nc
        ne   = Int(nEdges_per_cell[c])
        lmin = Inf; lmax = -Inf; lsum = 0.0
        for j in 1:ne
            l    = edge_lengths[cells_edges[c][j]]
            lmin = min(lmin, l)
            lmax = max(lmax, l)
            lsum += l
        end
        cell_d[c] = (lmax - lmin) / (lsum / ne)
    end
    return mean(cell_d)
end

function print_metrics(mesh, level)
    d  = distortion_metric(mesh)
    nc = length(mesh.cells.position)
    println("  Metrics (level $level): nc = $nc, mean cell irregularity = $(round(d, digits=4))")
end

function save_mesh_png(filename, mesh, label)
    fig = Figure(size = (700, 700))
    ax  = Axis(fig[1, 1], title = label, aspect = DataAspect())
    plotdualmesh!(ax, mesh)  # triangulation (orange) underneath
    plotmesh!(ax, mesh)      # Voronoi (blue) on top
    hidedecorations!(ax)
    GLMakie.save(filename, fig)
    return nothing
end

function save_level(outdir, mesh, level)
    nc    = length(mesh.cells.position)
    label = "mesh_periodic_regular_L$(level)_nc$(nc)"
    save_voronoi_to_vtu(joinpath(outdir, "$(label)_vor.vtu"), mesh)
    save_triangulation_to_vtu(joinpath(outdir, "$(label)_tri.vtu"), mesh)
    save_mesh_png(joinpath(outdir, "$(label).png"), mesh, "L$level — $nc cells")
    println("  Saved: $label")
    print_metrics(mesh, level)
end

function main(nc_ref, num_levels)
    outdir = "output"
    mkpath(outdir)

    dc = sqrt(X_PERIOD * Y_PERIOD / nc_ref)
    println("Level 0: creating regular hex mesh (reference nc≈$nc_ref, dc=$(round(dc, digits=4)))...")
    mesh = create_planar_hex_mesh(X_PERIOD, Y_PERIOD, dc)
    save_level(outdir, mesh, 0)

    for i in 1:num_levels
        nc_prev = length(mesh.cells.position)
        println("Level $i: refining from $nc_prev cells (cells + edge midpoints as generators)...")
        generators = vcat(mesh.cells.position, mesh.edges.position)
        mesh = VoronoiMesh(generators, X_PERIOD, Y_PERIOD)
        save_level(outdir, mesh, i)
    end

    println("\nDone! Output written to $(abspath(outdir))/")
    return nothing
end

nc_ref     = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 16
num_levels = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 4

main(nc_ref, num_levels)
