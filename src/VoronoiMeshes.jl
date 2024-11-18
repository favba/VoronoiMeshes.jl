module VoronoiMeshes

using Zeros, TensorsLite, ImmutableVectors, TensorsLiteGeometry

export @parallel

export VoronoiDiagram, PlanarVoronoiDiagram, SphericalVoronoiDiagram
export CellInfo, Cells, VertexInfo, Vertices, EdgeInfo, Edges
export VoronoiMesh, on_sphere, max_edges, integer_type, float_type, get_diagram
export meshplot, meshplot!, diagramplot, diagramplot!
export graph_partition, find_obtuse_triangles, periodic_edges_mask, periodic_vertices_mask
export check_mesh, check_edge_normal_and_tangent, check_vertex_indexing, check_cell_indexing

const VecMaybe1DxArray{TX, TYZ, N} = TensorsLite.VecArray{Vec{Union{TX, TYZ}, 1, TX, TYZ, TYZ}, N, Array{TX, N}, Array{TYZ, N}, Array{TYZ, N}}
const Vec1DxOr2DxyArray{TX, TXY, N} = TensorsLite.VecArray{Vec{Union{TX, Zero}, 1, TX, TXY, Zero}, N, Array{TX, N}, Array{TXY, N}, Array{Zero, N}}
const VecMaybe1DyArray{TY, TXZ, N} = TensorsLite.VecArray{Vec{Union{TY, TXZ}, 1, TXZ, TY, TXZ}, N, Array{TXZ, N}, Array{TY, N}, Array{TXZ, N}}

include("utils_pre.jl")

include("voronoi_diagrams.jl")

include("cell_struct.jl")

include("vertex_struct.jl")

include("edge_struct.jl")

struct VoronoiMesh{S, maxEdges, TI, TF, Tz}
    cells::Cells{S, maxEdges, TI, TF, Tz}
    vertices::Vertices{S, maxEdges, TI, TF, Tz}
    edges::Edges{S, maxEdges, TI, TF, Tz}

    function VoronoiMesh(
                c::Cells{S, NE, TI, TF, Tz},
                v::Vertices{S, NE, TI, TF, Tz},
                e::Edges{S, NE, TI, TF, Tz}
            ) where {S, NE, TI, TF, Tz}

        if !(get_diagram(c) === get_diagram(v) === get_diagram(e))
            throw(DimensionMismatch("`Cell`, `Vertices` and `Edges` structs are not based on the same Voronoi diagram"))
        end
        return new{S, NE, TI, TF, Tz}(c, v, e)
    end
end

get_diagram(mesh::VoronoiMesh) = get_diagram(mesh.cells)

Base.propertynames(::VoronoiMesh{false}) = (fieldnames(VoronoiMesh)..., :diagram, :x_period, :y_period)
Base.propertynames(::VoronoiMesh{true}) = (fieldnames(VoronoiMesh)..., :diagram, :sphere_radius)

Base.getproperty(mesh::VoronoiMesh, s::Symbol) = _getproperty(mesh, Val(s))
_getproperty(mesh::VoronoiMesh, ::Val{s}) where {s} = getfield(mesh, s)
_getproperty(mesh::VoronoiMesh, ::Val{:diagram}) = get_diagram(mesh)
_getproperty(mesh::VoronoiMesh{false}, ::Val{:x_period}) = get_diagram(mesh).x_period
_getproperty(mesh::VoronoiMesh{false}, ::Val{:y_period}) = get_diagram(mesh).y_period
_getproperty(mesh::VoronoiMesh{true}, ::Val{:sphere_radius}) = get_diagram(mesh).sphere_radius

for T in (
        :VoronoiMesh,
        :Cells, :CellInfo,
        :Vertices, :VertexInfo,
        :Edges, :EdgeInfo,
        :VoronoiDiagram,
    )
    @eval begin
        on_sphere(::Type{<:$T{S}}) where {S} = S
        max_edges(::Type{<:$T{S, N}}) where {S, N} = N
        integer_type(::Type{<:$T{S, N, TI}}) where {S, N, TI} = TI
        float_type(::Type{<:$T{S, N, TI, TF}}) where {S, N, TI, TF} = TF
    end
    for func in (:on_sphere, :max_edges, :integer_type, :float_type)
        @eval $func(field::$T) = $func(typeof(field))
    end
end

function VoronoiMesh(voronoi::VoronoiDiagram)
    verticesOnCell = voronoi.verticesOnCell
    cellsOnVertex = voronoi.cellsOnVertex
    cpos = voronoi.generators
    vpos = voronoi.vertices

    edges, vertex_pair_to_edge, cell_pair_to_edge = build_edges(voronoi)

    cellsOnCell = compute_cellsOnCell(verticesOnCell, cellsOnVertex)
    edgesOnCell = compute_edgesOnCell(cellsOnCell, cell_pair_to_edge)

    cells = Cells(length(cpos), cpos, verticesOnCell.length, verticesOnCell, edgesOnCell, cellsOnCell, CellInfo(voronoi))

    edgesOnVertex = compute_edgesOnVertex(cellsOnVertex, cell_pair_to_edge)

    vertices = Vertices(length(vpos), vpos, edgesOnVertex, cellsOnVertex, VertexInfo(voronoi))

    return VoronoiMesh(cells, vertices, edges)
end

include("utils_pos.jl")

#functions to be used when Makie is loaded
#Definitions are in ext/GeometryBasicsExt.jl
function create_cells_polygons_periodic end
function create_cells_polygons end
function create_dual_triangles_periodic end
function create_dual_triangles end
function create_edge_quadrilaterals_periodic end
function create_edge_quadrilaterals end
function create_cell_linesegments_periodic end
function create_cell_linesegments end
function create_diagram_linesegments_periodic end
function create_diagram_linesegments end

# Defined in ext/MakieExt.jl
function meshplot end
function meshplot! end

function diagramplot end
function diagramplot! end

end # Module
