# mesh_metrics.jl
#
# Per-cell quality metrics for periodic planar Voronoi meshes.
# Each metric function takes a VoronoiMesh and returns a Vector with one
# value per cell. Add new metrics here and register them in `METRICS` so
# they automatically show up in scripts that iterate over it.

module MeshMetrics

using TensorsLiteGeometry: closest
using Statistics: mean

export METRICS, cell_area, cell_area_normalized, cell_distortion, cell_distortion_rms,
       cell_diameter, cell_diameter_normalized, cell_alignment

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

end # module
