module VoronoiMeshes

using Zeros, TensorsLite, ImmutableVectors, TensorsLiteGeometry

export @parallel

export VoronoiDiagram, PlanarVoronoiDiagram, SphericalVoronoiDiagram
export CellInfo, Cells

include("utils_pre.jl")

include("voronoi_diagrams.jl")

include("cell_struct.jl")

struct VoronoiMesh{S, maxEdges, TI, TF, Tz}
end

end
