abstract type AbstractVoronoiDiagram{S, nEdges, TI <: Integer, TF <: Real, Tz} end

struct PlanarVoronoiDiagram{nEdges, TI, TF} <: AbstractVoronoiDiagram{false, nEdges, TI, TF, Zeros.Zero}
    generators::Vec2DxyArray{TF, 1}
    vertices::Vec2DxyArray{TF, 1}
    verticesOnCell::ImVecArray{nEdges, TI, 1}
    cellsOnVertex::Vector{NTuple{3, TI}}
    meshDensity::Vector{TF}
    x_period::TF
    y_period::TF
end

struct SphericalVoronoiDiagram{nEdges, TI, TF} <: AbstractVoronoiDiagram{true, nEdges, TI, TF, TF}
    generators::Vec3DArray{TF, 1}
    vertices::Vec3DArray{TF, 1}
    verticesOnCell::ImVecArray{nEdges, TI, 1}
    cellsOnVertex::Vector{NTuple{3, TI}}
    meshDensity::Vector{TF}
    sphere_radius::TF
end

struct VoronoiDiagram{S, nEdges, TI, TF, Tz} <: AbstractVoronoiDiagram{S, nEdges, TI, TF, Tz}
    data::Union{PlanarVoronoiDiagram{nEdges, TI, TF}, SphericalVoronoiDiagram{nEdges, TI, TF}}

    function VoronoiDiagram(d::PlanarVoronoiDiagram{nEdges, TI, TF}) where {nEdges, TI, TF}
        return new{false, nEdges, TI, TF, Zero}(d)
    end

    function VoronoiDiagram(d::SphericalVoronoiDiagram{nEdges, TI, TF}) where {nEdges, TI, TF}
        return new{true, nEdges, TI, TF, TF}(d)
    end
end

get_diagram(v::VoronoiDiagram{true, nEdges, TI, TF, TF}) where {nEdges, TI, TF} = Base.getfield(v,:data)::SphericalVoronoiDiagram{nEdges, TI, TF}
get_diagram(v::VoronoiDiagram{false, nEdges, TI, TF, Zero}) where {nEdges, TI, TF} = Base.getfield(v,:data)::PlanarVoronoiDiagram{nEdges, TI, TF}

Base.getproperty(v::VoronoiDiagram, s::Symbol) = _getproperty(v, Val{s}())
_getproperty(v::VoronoiDiagram, ::Val{:data}) = get_diagram(v)

for s in fieldnames(PlanarVoronoiDiagram)
    @eval _getproperty(v::VoronoiDiagram{false}, ::Val{$(QuoteNode(s))}) = getfield(get_diagram(v), $(QuoteNode(s)))
end

for s in fieldnames(SphericalVoronoiDiagram)
    @eval _getproperty(v::VoronoiDiagram{true}, ::Val{$(QuoteNode(s))}) = getfield(get_diagram(v), $(QuoteNode(s)))
end
