mutable struct VertexInfo{S, NE, TI, TF, Tz}
    const diagram::VoronoiDiagram{S, NE, TI, TF, Tz}
    centroid::TensorsLite.VecMaybe2DxyArray{TF, Tz, 1}
    area::Vector{TF}
    kiteAreas::Vector{NTuple{3, TF}}
    longitude::Vector{TF}
    latitude::Vector{TF}

    function VertexInfo(diagram::VoronoiDiagram{S, NE, TI, TF, Tz}) where {S, NE, TI, TF, Tz}
        return new{S, NE, TI, TF, Tz}(diagram)
    end
end

struct Vertices{S, NE, TI, TF, Tz}
    n::Int
    position::TensorsLite.VecMaybe2DxyArray{TF,Tz, 1}
    edges::Vector{NTuple{3, TI}}
    cells::Vector{NTuple{3, TI}}
    info::VertexInfo{S, NE, TI, TF, Tz}
end

Base.getproperty(cell::Vertices, s::Symbol) = _getproperty(cell, Val(s))
_getproperty(cell::Vertices, ::Val{s}) where s = getfield(cell, s)
_getproperty(cell::Vertices{false}, ::Val{:x_period}) = getfield(cell, :info).diagram.x_period
_getproperty(cell::Vertices{false}, ::Val{:y_period}) = getfield(cell, :info).diagram.y_period
_getproperty(cell::Vertices{true}, ::Val{:sphere_radius}) = getfield(cell, :info).diagram.sphere_radius

for s in fieldnames(VertexInfo)
    func = Symbol(string("compute_vertex_",s))
    @eval function _getproperty(vertex::Vertices, ::Val{$(QuoteNode(s))})
        info = getfield(vertex, :info)
        if !isdefined(info, $(QuoteNode(s)))
            setfield!(info, $(QuoteNode(s)), $func(vertex))
        end
        return getfield(info, $(QuoteNode(s)))
    end
end

include("vertex_info_creation.jl")

