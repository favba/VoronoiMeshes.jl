
mutable struct CellInfo{S, nEdges, TI, TF, Tz}
    const diagram::VoronoiDiagram{S, nEdges, TI, TF, Tz}
    centroid::TensorsLite.VecMaybe2DxyArray{TF, Tz, 1}
    area::Vector{TF}
    longitude::Vector{TF}
    latitude::Vector{TF}
    verticalUnitVectors::TensorsLite.VecMaybe2DxyArray{Tz, TF, 1} # Actually, Maybe1DzArray
    tangentPlane::NTuple{2, TensorsLite.VecMaybe2DxyArray{TF, Tz, 1}}

    function CellInfo(diagram::VoronoiDiagram{S, nEdges, TI, TF, Tz}) where {S, nEdges, TI, TF, Tz}
        return new{S, nEdges, TI, TF, Tz}(diagram)
    end
end

struct Cells{S, nEdges, TI, TF, Tz}
    n::Int
    position::TensorsLite.VecMaybe2DxyArray{TF,Tz, 1}
    nEdges::Vector{UInt8}
    vertices::ImVecArray{nEdges, TI, 1}
    edges::ImVecArray{nEdges, TI, 1}
    cells::ImVecArray{nEdges, TI, 1}
    info::CellInfo{S, nEdges, TI, TF, Tz}
end

Base.getproperty(cell::Cells, s::Symbol) = _getproperty(cell, Val(s))
_getproperty(cell::Cells, ::Val{s}) where s = getfield(cell, s)
_getproperty(cell::Cells{false}, ::Val{:x_period}) = getfield(cell, :info).diagram.x_period
_getproperty(cell::Cells{false}, ::Val{:y_period}) = getfield(cell, :info).diagram.y_period
_getproperty(cell::Cells{true}, ::Val{:sphere_radius}) = getfield(cell, :info).diagram.sphere_radius

for s in fieldnames(CellInfo)
    func = Symbol(string("compute_cell_",s))
    @eval function _getproperty(cell::Cells, ::Val{$(QuoteNode(s))})
        info = getfield(cell, :info)
        if !isdefined(info, $(QuoteNode(s)))
            setfield!(info, $(QuoteNode(s)), $func(cell))
        end
        return getfield(info, $(QuoteNode(s)))
    end
end

include("cell_info_creation.jl")
