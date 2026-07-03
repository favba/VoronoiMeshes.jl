# build_set_irregular_meshes.jl
#
# Generates a series of nc-cell Voronoi meshes with increasing irregularity.
#   Level 0 : fully converged centroidal mesh (most regular)
#   Levels 1..num_levels : each level takes the previous level's generator
#              points, adds growing random noise, then runs PERTURB_ITERS
#              Lloyd steps — enough to fix degenerate cells without erasing
#              the introduced irregularity.
#
# Usage:
#   julia --project=. build_set_irregular_meshes.jl [nc] [num_levels] [base_strength]
#
# Defaults: nc=80, num_levels=5, base_strength=0.15
#
# base_strength is the perturbation amplitude at level 1, expressed as a
# fraction of the average cell spacing dc ≈ 1/√nc.  It scales linearly:
# level i uses strength = base_strength * i.
#
# Note: PNG export uses GLMakie. On headless servers swap to CairoMakie
# (add it to the environment and replace `using GLMakie` below).

using VoronoiMeshes
using DelaunayTriangulation
using TensorsLite
using WriteVTK
using GLMakie
using Statistics

const PERTURB_ITERS = 5
const X_PERIOD = 1.0
const Y_PERIOD = 1.0

# Perturbation band: only cells within BAND_WIDTH of the line y = LINE_SLOPE*x + LINE_INTERCEPT
# are perturbed. Cells outside the band keep their positions unchanged.
const LINE_SLOPE = 0.3
const LINE_INTERCEPT = 0.3
const BAND_WIDTH = 0.12

# Extra Lloyd iterations applied one at a time after perturbation to eliminate
# obtuse Delaunay triangles (circumcenter outside triangle). Set to 0 to skip.
const MAX_FIX_ITERS = 10

# Per-cell irregularity: (max_edge - min_edge) / mean_edge for each cell,
# then averaged globally. Zero for perfectly regular cells, grows with distortion.
function distortion_metric(mesh)
    edge_lengths = mesh.edges.length
    cells_edges = mesh.cells.edges
    nEdges_per_cell = mesh.cells.nEdges
    nc = length(nEdges_per_cell)

    cell_d = Vector{Float64}(undef, nc)
    for c in 1:nc
        ne = Int(nEdges_per_cell[c])
        lmin = Inf;
        lmax = -Inf;
        lsum = 0.0
        for j in 1:ne
            l = edge_lengths[cells_edges[c][j]]
            lmin = min(lmin, l)
            lmax = max(lmax, l)
            lsum += l
        end
        cell_d[c] = (lmax - lmin) / (lsum / ne)
    end
    return mean(cell_d)
end

function print_metrics(mesh, level)
    d = distortion_metric(mesh)
    n_obtuse = length(find_obtuse_triangles(mesh))
    n_tri = length(mesh.vertices.cells)
    println("  Distortion (level $level): mean cell irregularity = $(round(d, digits=4)), obtuse triangles = $n_obtuse / $n_tri, x_period = $(mesh.x_period), y_period = $(mesh.y_period)")
end

function save_mesh_png(filename, mesh, label)
    fig = Figure(size=(700, 700))
    ax = Axis(fig[1, 1], title=label, aspect=DataAspect())
    plotdualmesh!(ax, mesh)  # triangulation (orange) underneath
    plotmesh!(ax, mesh)      # Voronoi (blue) on top
    hidedecorations!(ax)
    GLMakie.save(filename, fig)
    return nothing
end

function perturb_positions(mesh, strength)
    pos = mesh.cells.position
    nc = length(pos)
    dc = sqrt(X_PERIOD * Y_PERIOD / nc)

    # Perpendicular distance from (x,y) to line y = LINE_SLOPE*x + LINE_INTERCEPT,
    # i.e. -LINE_SLOPE*x + y - LINE_INTERCEPT = 0
    norm_factor = sqrt(1 + LINE_SLOPE^2)

    x_new = copy(pos.x)
    y_new = copy(pos.y)
    for c in 1:nc
        p = pos[c]
        dist = abs(p.y - LINE_SLOPE * p.x - LINE_INTERCEPT) / norm_factor
        if dist < BAND_WIDTH
            x_new[c] = mod(p.x + strength * dc * randn(), X_PERIOD)
            y_new[c] = mod(p.y + strength * dc * randn(), Y_PERIOD)
        end
    end
    return VecArray(x=x_new, y=y_new)
end

function main(nc, num_levels, base_strength)
    outdir = "output"
    mkpath(outdir)

    dc = sqrt(X_PERIOD * Y_PERIOD / nc)
    println("Level 0: creating regular hex mesh (reference nc≈$nc, dc=$(round(dc, digits=4)))...")
    hex_mesh = create_planar_hex_mesh(X_PERIOD, Y_PERIOD, dc)
    # create_planar_hex_mesh rounds the cell counts to fit an integer number of hex
    # rows/columns, so its returned domain isn't exactly X_PERIOD x Y_PERIOD. Reuse its
    # generators (same count) but rebuild against the exact periodic domain, which lets
    # Lloyd relaxation (the default) spread them to fill it.
    mesh = VoronoiMesh(hex_mesh.cells.position, X_PERIOD, Y_PERIOD)
    actual_nc = length(mesh.cells.position)
    println("  Actual cell count: $actual_nc")

    label = "mesh_periodic_irregular_L0_nc$(actual_nc)_d0.0"
    save_voronoi_to_vtu(joinpath(outdir, "$(label)_vor.vtu"), mesh)
    save_triangulation_to_vtu(joinpath(outdir, "$(label)_tri.vtu"), mesh)
    save_mesh_png(joinpath(outdir, "$(label).png"), mesh, label)
    println("  Saved: $label")
    print_metrics(mesh, 0)

    prev_mesh = mesh
    for i in 1:num_levels
        strength = base_strength * i
        println("Level $i: perturbing (strength=$(round(strength, digits=3))), $PERTURB_ITERS Lloyd iterations...")
        generators = perturb_positions(prev_mesh, strength)
        mesh = VoronoiMesh(generators, X_PERIOD, Y_PERIOD; max_iter=PERTURB_ITERS)

        # Fix any obtuse triangles with extra Lloyd iterations (one at a time)
        extra = 0
        while !isempty(find_obtuse_triangles(mesh)) && extra < MAX_FIX_ITERS
            mesh = VoronoiMesh(mesh.cells.position, X_PERIOD, Y_PERIOD; max_iter=4)
            extra += 1
        end
        extra > 0 && println("  Applied $extra extra Lloyd iteration(s) to remove obtuse triangles.")

        label = "mesh_periodic_irregular_L$(i)_nc$(actual_nc)_d$(round(strength, digits=3))"
        save_voronoi_to_vtu(joinpath(outdir, "$(label)_vor.vtu"), mesh)
        save_triangulation_to_vtu(joinpath(outdir, "$(label)_tri.vtu"), mesh)
        save_mesh_png(joinpath(outdir, "$(label).png"), mesh, label)
        println("  Saved: $label")
        print_metrics(mesh, i)

        prev_mesh = mesh
    end

    println("\nDone! Output written to $(abspath(outdir))/")
    return nothing
end

nc = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 80
num_levels = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 5
base_strength = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 0.15

main(nc, num_levels, base_strength)
