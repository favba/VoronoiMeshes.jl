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
    position::TensorsLite.VecMaybe2DxyArray{TF, Tz, 1}
    edges::Vector{NTuple{3, TI}}
    cells::Vector{NTuple{3, TI}}
    info::VertexInfo{S, NE, TI, TF, Tz}
end

Base.getproperty(cell::Vertices, s::Symbol) = _getproperty(cell, Val(s))
_getproperty(cell::Vertices, ::Val{s}) where {s} = getfield(cell, s)
_getproperty(cell::Vertices{false}, ::Val{:x_period}) = getfield(cell, :info).diagram.x_period
_getproperty(cell::Vertices{false}, ::Val{:y_period}) = getfield(cell, :info).diagram.y_period
_getproperty(cell::Vertices{true}, ::Val{:sphere_radius}) = getfield(cell, :info).diagram.sphere_radius

include("vertex_info_creation.jl")

for s in fieldnames(VertexInfo)
    if s !== :diagram
        func = Symbol(string("compute_vertex_", s))
        @eval function _getproperty(vertex::Vertices, ::Val{$(QuoteNode(s))})
            info = getfield(vertex, :info)
            if !isdefined(info, $(QuoteNode(s)))
                setfield!(info, $(QuoteNode(s)), $func(vertex))
            end
            return getfield(info, $(QuoteNode(s)))
        end
        for nEdges in 6:10
            for TI in (Int32, Int64)
                for TF in (Float32, Float64)
                    @eval precompile($func, (Vertices{true, $nEdges, $TI, $TF, $TF},))
                    @eval precompile($func, (Vertices{false, $nEdges, $TI, $TF, Zero},))
                end
            end
        end
    else
        _getproperty(vertex::Vertices, ::Val{:diagram}) = getfield(vertex, :info).diagram
    end
end
