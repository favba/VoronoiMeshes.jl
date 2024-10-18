module VoronoiMeshes

using Zeros, TensorsLite, ImmutableVectors, TensorsLiteGeometry

export @parallel

export VoronoiDiagram, PlanarVoronoiDiagram, SphericalVoronoiDiagram

include("utils_pre.jl")

include("voronoi_diagrams.jl")

end
