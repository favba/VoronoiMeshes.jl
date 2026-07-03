# mesh_metrics.jl
#
# Per-cell quality metrics for periodic planar Voronoi meshes.
# Each metric function takes a VoronoiMesh and returns a Vector with one
# value per cell. Add new metrics here and register them in `METRICS` so
# they automatically show up in scripts that iterate over it.

module MeshMetrics

using TensorsLiteGeometry: closest

export METRICS, cell_area, cell_distortion, cell_diameter

# Cell area, already computed and cached by VoronoiMeshes.
cell_area(mesh) = mesh.cells.area

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

# Ordered (name, function) pairs used by scripts to compute and report
# every available per-cell metric without hardcoding the list.
const METRICS = (
    ("area", cell_area),
    ("distortion", cell_distortion),
    ("diameter", cell_diameter),
)

end # module
