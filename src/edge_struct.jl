mutable struct EdgeInfo{S, NE, TI, TF, Tz}
    const diagram::VoronoiDiagram{S, NE, TI, TF, Tz}
    midpoint::TensorsLite.VecMaybe2DxyArray{TF, Tz, 1}
    length::Vector{TF}
    cellsDistance::Vector{TF}
    angle::Vector{TF}
    longitude::Vector{TF}
    latitude::Vector{TF}
    normal::TensorsLite.VecMaybe2DxyArray{TF, Tz, 1}
    tangent::TensorsLite.VecMaybe2DxyArray{TF, Tz, 1}

    function EdgeInfo(diagram::VoronoiDiagram{S, NE, TI, TF, Tz}) where {S, NE, TI, TF, Tz}
        return new{S, NE, TI, TF, Tz}(diagram)
    end
end

struct Edges{S, NE, TI, TF, Tz}
    n::Int
    position::TensorsLite.VecMaybe2DxyArray{TF,Tz, 1}
    vertices::Vector{NTuple{2, TI}}
    cells::Vector{NTuple{2, TI}}
    info::EdgeInfo{S, NE, TI, TF, Tz}
end

Base.getproperty(edge::Edges, s::Symbol) = _getproperty(edge, Val(s))
_getproperty(edge::Edges, ::Val{s}) where s = getfield(edge, s)
_getproperty(edge::Edges{false}, ::Val{:x_period}) = getfield(edge, :info).diagram.x_period
_getproperty(edge::Edges{false}, ::Val{:y_period}) = getfield(edge, :info).diagram.y_period
_getproperty(edge::Edges{true}, ::Val{:sphere_radius}) = getfield(edge, :info).diagram.sphere_radius

for s in fieldnames(EdgeInfo)
    func = Symbol(string("compute_edge_",s))
    @eval function _getproperty(edge::Edges, ::Val{$(QuoteNode(s))})
        info = getfield(edge, :info)
        if !isdefined(info, $(QuoteNode(s)))
            setfield!(info, $(QuoteNode(s)), $func(edge))
        end
        return getfield(info, $(QuoteNode(s)))
    end
end

include("edge_info_creation.jl")

function Edges(diagram::VoronoiDiagram{false, NE, TI, TF}) where {NE, TI, TF}

    cpos = diagram.generators
    verticesOnCell = diagram.verticesOnCell
    cellsOnVertex = diagram.cellsOnVertex
    nVertex = length(cellsOnVertex)
    nCells = length(verticesOnCell)
    nEdges = nVertex + nCells

    position = similar(diagram.vertices, nEdges)
    vertices = Vector{NTuple{2, TI}}(undef, nEdges)
    cells = Vector{NTuple{2, TI}}(undef, nEdges)

    xp = diagram.x_period
    yp = diagram.y_period

    touched_vertex_pair = Set{NTuple{2,TI}}()
    e = 0
    @inbounds for c in eachindex(cpos)
        cell_vertices = verticesOnCell[c]
        cp = cpos[c]
        v_1 = cell_vertices[1]
        v2 = v_1
        for v in 2:length(cell_vertices)
            v1 = v2
            v2 = cell_vertices[v]
            pair = ordered(v1,v2)
            if !(pair in touched_vertex_pair)
                e+=one(TI)
                push!(touched_vertex_pair, pair)
                c2 = find_cellOnCell(c, v2, v1, cellsOnVertex)
                position[e] = periodic_to_base_point((cp + closest(cp, cpos[c2], xp, yp) / 2), xp, yp)
                vertices[e] =  (v1, v2)
                cells[e] = (c, c2)
            end
        end
        v1 = v2
        v2 = v_1
        pair = ordered(v1,v2)
        if !(pair in touched_vertex_pair)
            e+=one(TI)
            push!(touched_vertex_pair, pair)
            c2 = find_cellOnCell(c, v2, v1, cellsOnVertex)
            position[e] = periodic_to_base_point((cp + closest(cp, cpos[c2], xp, yp) / 2), xp, yp)
            vertices[e] =  (v1, v2)
            cells[e] = (c, c2)
        end
    end

    @assert e == nEdges

    return Edges(e, position, vertices, cells, EdgeInfo(diagram))
end

function Edges(diagram::VoronoiDiagram{true, NE, TI, TF}) where {NE, TI, TF}

    cpos = diagram.generators
    verticesOnCell = diagram.verticesOnCell
    cellsOnVertex = diagram.cellsOnVertex
    nVertex = length(cellsOnVertex)
    nCells = length(verticesOnCell)
    nEdges = nVertex + nCells - 2

    position = similar(diagram.vertices, nEdges)
    vertices = Vector{NTuple{2, TI}}(undef, nEdges)
    cells = Vector{NTuple{2, TI}}(undef, nEdges)

    R = diagram.sphere_radius

    touched_vertex_pair = Set{NTuple{2,TI}}()
    e = 0
    @inbounds for c in eachindex(cpos)
        cell_vertices = verticesOnCell[c]
        cp = cpos[c]
        v_1 = cell_vertices[1]
        v2 = v_1
        for v in 2:length(cell_vertices)
            v1 = v2
            v2 = cell_vertices[v]
            pair = ordered(v1,v2)
            if !(pair in touched_vertex_pair)
                e+=one(TI)
                push!(touched_vertex_pair, pair)
                c2 = find_cellOnCell(c, v2, v1, cellsOnVertex)
                position[e] = arc_midpoint(R, cp, cpos[c2])
                vertices[e] =  (v1, v2)
                cells[e] = (c, c2)
            end
        end
        v1 = v2
        v2 = v_1
        pair = ordered(v1,v2)
        if !(pair in touched_vertex_pair)
            e+=one(TI)
            push!(touched_vertex_pair, pair)
            c2 = find_cellOnCell(c, v2, v1, cellsOnVertex)
            position[e] = arc_midpoint(R, cp, cpos[c2])
            vertices[e] =  (v1, v2)
            cells[e] = (c, c2)
        end
    end

    @assert e == nEdges

    return Edges(e, position, vertices, cells, EdgeInfo(diagram))
end
