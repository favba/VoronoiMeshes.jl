mutable struct CellInfo{S, max_nEdges, TI, TF, Tz}
    const diagram::VoronoiDiagram{S, max_nEdges, TI, TF, Tz}
    centroid::TensorsLite.VecMaybe2DxyArray{TF, Tz, 1}
    area::Vector{TF}
    longitude::Vector{TF}
    latitude::Vector{TF}
    normal::TensorsLite.VecMaybe2DxyArray{Tz, TF, 1} # The same as Maybe1DzArray{TF, Tz, 1}
    zonalVector::Vec1DxOr2DxyArray{TF, Tz, 1}
    meridionalVector::VecMaybe1DyArray{TF, Tz, 1}

    function CellInfo(diagram::VoronoiDiagram{S, max_nEdges, TI, TF, Tz}) where {S, max_nEdges, TI, TF, Tz}
        return new{S, max_nEdges, TI, TF, Tz}(diagram)
    end
end

const planar_cellinfo_names = (:centroid, :area, :x_period, :y_period)
const spherical_cellinfo_names = (filter(!=(:diagram), fieldnames(CellInfo))..., :sphere_radius)

for nEdges in 6:10
    precompile(CellInfo, (VoronoiDiagram{false, nEdges, Int32, Float64, Zero},))
    precompile(CellInfo, (VoronoiDiagram{true, nEdges, Int32, Float64, Float64},))
end

struct Cells{S, max_nEdges, TI, TF, Tz}
    n::Int
    position::TensorsLite.VecMaybe2DxyArray{TF, Tz, 1}
    nEdges::Vector{UInt8}
    vertices::ImVecArray{max_nEdges, TI, 1}
    edges::ImVecArray{max_nEdges, TI, 1}
    cells::ImVecArray{max_nEdges, TI, 1}
    info::CellInfo{S, max_nEdges, TI, TF, Tz}
end

get_diagram(c::Cells) = getfield(c, :info).diagram

for nEdges in 6:10
    precompile(
        Cells, (
            Int, Vec2DxyArray{Float64, 1}, Vector{UInt8}, ImVecArray{nEdges, Int32, 1},
            ImVecArray{nEdges, Int32, 1}, ImVecArray{nEdges, Int32, 1},
            CellInfo{false, nEdges, Int32, Float64, Zero},
        )
    )
    precompile(
        Cells, (
            Int, Vec3DArray{Float64, 1}, Vector{UInt8}, ImVecArray{nEdges, Int32, 1},
            ImVecArray{nEdges, Int32, 1}, ImVecArray{nEdges, Int32, 1},
            CellInfo{true, nEdges, Int32, Float64, Float64},
        )
    )
end

const cell_names = (:n, :position, :nEdges, :vertices, :edges, :cells)

Base.propertynames(::Cells{false}) = (cell_names..., planar_cellinfo_names...)
Base.propertynames(::Cells{true}) = (cell_names..., spherical_cellinfo_names...)

Base.getproperty(cell::Cells, s::Symbol) = _getproperty(cell, Val(s))
_getproperty(cell::Cells, ::Val{s}) where {s} = getfield(cell, s)
_getproperty(cell::Cells{false}, ::Val{:x_period}) = get_diagram(cell).x_period
_getproperty(cell::Cells{false}, ::Val{:y_period}) = get_diagram(cell).y_period
_getproperty(cell::Cells{true}, ::Val{:sphere_radius}) = get_diagram(cell).sphere_radius

include("cell_info_creation.jl")

for s in fieldnames(CellInfo)
    if s !== :diagram
        func = Symbol(string("compute_cell_", s))
        @eval function _getproperty(cell::Cells, ::Val{$(QuoteNode(s))})
            info = getfield(cell, :info)
            if !isdefined(info, $(QuoteNode(s)))
                setfield!(info, $(QuoteNode(s)), $func(cell))
            end
            return getfield(info, $(QuoteNode(s)))
        end

        for nEdges in 6:10
            for TI in (Int32, Int64)
                for TF in (Float32, Float64)
                    @eval precompile($func, (Cells{true, $nEdges, $TI, $TF, $TF},))
                    @eval precompile($func, (Cells{false, $nEdges, $TI, $TF, Zero},))
                end
            end
        end
    else
        _getproperty(cell::Cells, ::Val{:diagram}) = getfield(cell, :info).diagram
    end
end
