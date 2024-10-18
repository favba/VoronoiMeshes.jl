module VoronoiMeshes

using TensorsLite, ImmutableVectors, TensorsLiteGeometry

export @parallel

export VoronoiDiagram, PlanarPeriodicVoronoiDiagram, SphericalVoronoiDiagram

include("utils_pre.jl")

include("voronoi_diagrams.jl")

end
