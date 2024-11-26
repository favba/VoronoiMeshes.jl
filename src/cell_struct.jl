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

"""
    Cells{OnSphere, max_nEdges, <:Integer, <:Float, <:Union{Float, Zeros.Zero}} 

Struct that holds all Voronoi Cells information in a struct of arrays (SoA) layout, that is,
each field is usually an array with the requested data for each cell.

Data should be accessed using the `getproperty` function, preferably through the "dot" syntax.
For example, to fetch an array with each cells area, use `cells.area`.

When the struct is initialized only part of its data is actually existent. We refer to
those fields as the "Base Data". Other fields are only computed and stored if they are called at least once.
We refer to those as the "Computed Data".

# Base Data
- `n::Int`: The total number of Voronoi Cells.
- `position::VecArray`: An array with each cells position vector, that is, the position of the Voronoi
   Cell generator point. An array with a particular coordinate can also be extracted throught the dot
   syntax. For example, an array with `x` coordinates of the cells is given by `cells.position.x`.
- `nEdges::Vector{UInt8}`: The number of edges on each cell (which is equal to the number of vertices).
- `vertices::ImmutableVectorArray`: A Vector of vectors with the index ID of the vertices forming the cell.
- `edges::ImmutableVectorArray`: A Vector of vectors with the index ID of the edges forming the cell.
- `cells::ImmutableVectorArray`: A Vector of vectors with the index ID of the cells neighboring a given cell.
- `sphere_radius::Real` (Spherical meshes only): The sphere radius.
- `x_period::Real` (Planar meshes only): The domain `x` direction period.
- `y_period::Real` (Planar meshes only): The domain `y` direction period.

# Computed Data
- `area::Vector`: An array with the area of each Voronoi cell.
- `centroid::VecArray`: An array with each cells centroid position vector. For Centroidal Voronoi
  meshes with constant density function this should virtually be the same as the `position` vector.
- `longitude::Vector`(Spherical meshes only): The longitude in radians of the cell `position` vector.
- `latitude::Vector`(Spherical meshes only): The latitude in radians of the cell `position` vector.
- `normal::VecArray`(Spherical meshes only): The unit vector perpendicular to the plane tangent to the
  sphere at the cell `position`.
- `zonalVector::VecArray`(Spherical meshes only): The unit vector tangent to the sphere and pointing eastward.
- `meridionalVector::VecArray`(Spherical meshes only): The unit vector tangent to the sphere and pointing northward.

"""
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
