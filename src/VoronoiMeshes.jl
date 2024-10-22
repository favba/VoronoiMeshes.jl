module VoronoiMeshes

using Zeros, TensorsLite, ImmutableVectors, TensorsLiteGeometry

export @parallel

export VoronoiDiagram, PlanarVoronoiDiagram, SphericalVoronoiDiagram
export CellInfo, Cells, VertexInfo, Vertices, EdgeInfo, Edges
export VoronoiMesh

const VecMaybe1DxArray{TX, TYZ, N} = TensorsLite.VecArray{Union{TX, TYZ}, N, Array{TX, N}, Array{TYZ, N}, Array{TYZ, N}}
const VecMaybe1DyArray{TY, TXZ, N} = TensorsLite.VecArray{Union{TY, TXZ}, N, Array{TXZ, N}, Array{TY, N}, Array{TXZ, N}}

include("utils_pre.jl")

include("voronoi_diagrams.jl")

include("cell_struct.jl")

include("vertex_struct.jl")

include("edge_struct.jl")

struct VoronoiMesh{S, maxEdges, TI, TF, Tz}
    cells::Cells{S, maxEdges, TI, TF, Tz}
    vertices::Vertices{S, maxEdges, TI, TF, Tz}
    edges::Edges{S, maxEdges, TI, TF, Tz}
end

Base.getproperty(mesh::VoronoiMesh, s::Symbol) = _getproperty(mesh, Val(s))
_getproperty(mesh::VoronoiMesh, ::Val{s}) where s = getfield(mesh, s)
_getproperty(mesh::VoronoiMesh{false}, ::Val{:x_period}) = getfield(mesh, :cells).x_period
_getproperty(mesh::VoronoiMesh{false}, ::Val{:y_period}) = getfield(mesh, :cells).y_period
_getproperty(mesh::VoronoiMesh{true}, ::Val{:sphere_radius}) = getfield(mesh, :cells).sphere_radius

end
