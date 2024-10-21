module VoronoiMeshes

using Zeros, TensorsLite, ImmutableVectors, TensorsLiteGeometry

export @parallel

export VoronoiDiagram, PlanarVoronoiDiagram, SphericalVoronoiDiagram
export CellInfo, Cells, VertexInfo, Vertices

include("utils_pre.jl")

include("voronoi_diagrams.jl")

include("cell_struct.jl")

include("vertex_struct.jl")

struct VoronoiMesh{S, maxEdges, TI, TF, Tz}
    cells::Cells{S, maxEdges, TI, TF, Tz}
    vertices::Vertices{S, maxEdges, TI, TF, Tz}
end

end
