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

