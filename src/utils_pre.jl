"""
Works just like `Threads.@threads` (except you can't interpolate values with `\$`), but runs serially if `Threads.nthreads()==1` to avoid task creation overhead.
"""
macro parallel(ex)
    pex = quote
        let
            if $Threads.nthreads() == 1
                $(ex)
            else
                $Threads.@threads $(ex)
            end
        end
    end
    return esc(pex)
end 

function tmap!(output, func::F, var::Vararg) where F <: Function
    @parallel for i in eachindex(output)
        @inbounds output[i] = @inline func(map(x->@inbounds(x[i]), var)...)
    end
    return output
end

@inline function unsafe_drop_element(t::NTuple{3}, el)
    if t[1] == el
        return (t[2], t[3])
    elseif t[2] == el
        return (t[1], t[3])
    else
        return (t[1], t[2])
    end
end

@inline function unsafe_tuple2_intersection(t1::NTuple{2}, t2::NTuple{2})
    t11, t12 = t1
    t21, t22 = t2
    t11 == t21 ? t11 : t11 == t22 ? t11 : t12
end

@inline function find_cellOnCell(current_cell, shared_vertex1, shared_vertex2, cellsOnVertex)
    cells_v1 = cellsOnVertex[shared_vertex1]
    cells_v2 = cellsOnVertex[shared_vertex2]
    cells1 = unsafe_drop_element(cells_v1, current_cell)
    cells2 = unsafe_drop_element(cells_v2, current_cell)
    return unsafe_tuple2_intersection(cells1, cells2)
end

@inline ordered(x,y) = y < x ? (y, x) : (x, y)

function compute_cellsOnCell!(cellsOnCell, verticesOnCell, cellsOnVertex)
# The ordering of the indices is the same as seen in the meshes provided by NCAR
# They DO NOT follow what is described in their Mesh Specification document
# It seems the document (at least v1.0) specification is not correct.
# The ordering here is that cellsOncell[n] is between verticesOnCell[n+1] and verticesOnCell[n]
# Their document says cellsOnCell[n] should be between verticesOnCell[n] and verticesOnCell[n-1]
# but their own meshes don't follow this rule and uses the ordering used by the code below

    cellsOnCellTuple = cellsOnCell.data
    @parallel for c in eachindex(verticesOnCell)

        @inbounds begin
            cells = eltype(cellsOnCell)()
            vertices = verticesOnCell[c]
            lv = length(vertices)

            vlead = vertices[1]
            for i in 2:lv
                vprev = vlead
                vlead = vertices[i]
                cells = push(cells, find_cellOnCell(c, vlead, vprev, cellsOnVertex))
            end
            vprev = vlead
            vlead = vertices[1]
            cells = push(cells, find_cellOnCell(c, vlead, vprev, cellsOnVertex))

            cellsOnCellTuple[c] = cells.data
        end
    end
    return cellsOnCell
end

function compute_cellsOnCell(verticesOnCell, cellsOnVertex)
    cellsOnCell = ImmutableVectorArray(similar(verticesOnCell.data), verticesOnCell.length)
    return compute_cellsOnCell!(cellsOnCell, verticesOnCell, cellsOnVertex)
end

function compute_edgesOnCell!(edgesOnCell, cellsOnCell, cellsOnEdge::Vector{NTuple{2, TI}}) where TI

    edgesOnCellTuple = edgesOnCell.data
    @parallel for c in eachindex(cellsOnCell)
        @inbounds begin
            cells = cellsOnCell[c]
            edges = eltype(cellsOnCell)()
            for e in eachindex(cells)
                pair = ordered(TI(c), cells[e])
                e_i = findfirst(x -> (ordered(x[1], x[2]) === pair), cellsOnEdge)
                edges = push(edges, e_i)
            end
            edgesOnCellTuple[c] = edges.data
        end
    end
    return edgesOnCell
end

function compute_edgesOnCell(cellsOnCell, cellsOnEdge)
    edgesOnCell = ImmutableVectorArray(similar(cellsOnCell.data), cellsOnCell.length)
    return compute_edgesOnCell!(edgesOnCell, cellsOnCell, cellsOnEdge)
end

function compute_edgesOnVertex!(edgesOnVertex, cellsOnVertex, cellsOnEdge::Vector{NTuple{2, TI}}) where TI

    @parallel for v in eachindex(edgesOnVertex)
        @inbounds begin
            c1, c2, c3 = cellsOnVertex[v]
            pair1 = ordered(c3, c1)
            e1 = findfirst(x -> (ordered(x[1], x[2]) === pair1), cellsOnEdge)
            pair2 = ordered(c1, c2)
            e2 = findfirst(x -> (ordered(x[1], x[2]) === pair2), cellsOnEdge)
            pair3 = ordered(c2, c3)
            e3 = findfirst(x -> (ordered(x[1], x[2]) === pair3), cellsOnEdge)
            edgesOnVertex[v] = (e1, e2, e3)
       end
    end
    return edgesOnVertex
end

function compute_edgesOnVertex(cellsOnVertex::Vector{NTuple{3, TI}}, cellsOnEdge) where TI
    edgesOnVertex = Vector{NTuple{3,TI}}(undef,length(cellsOnVertex))
    return compute_edgesOnVertex!(edgesOnVertex, cellsOnVertex, cellsOnEdge)
end

function compute_polygon_area_periodic!(output,vpos,verticesOnPolygon,xp::Number,yp::Number)
    @parallel for c in eachindex(verticesOnPolygon)
        @inbounds output[c] = area(vpos,verticesOnPolygon[c],xp,yp)
    end
    return output
end

function compute_polygon_area_periodic(vpos,verticesOnPolygon,xp::Number,yp::Number)
    output = Vector{nonzero_eltype(eltype(vpos))}(undef,length(verticesOnPolygon))
    return compute_polygon_area_periodic!(output,vpos,verticesOnPolygon,xp,yp)
end

function compute_polygon_area_spherical!(output, vpos, verticesOnPolygon, R::Number)
    @parallel for c in eachindex(verticesOnPolygon)
        @inbounds output[c] = spherical_polygon_area(R, vpos,verticesOnPolygon[c])
    end
    return output
end

function compute_polygon_area_spherical(vpos, verticesOnPolygon, R::Number)
    output = Vector{nonzero_eltype(eltype(vpos))}(undef, length(verticesOnPolygon))
    return compute_polygon_area_spherical!(output, vpos, verticesOnPolygon, R)
end

function compute_polygon_centroid_periodic!(output,vpos,verticesOnPolygon,xp::Number,yp::Number)
    @parallel for c in eachindex(verticesOnPolygon)
        @inbounds output[c] = centroid(vpos,verticesOnPolygon[c],xp,yp)
    end
    return output
end

function compute_polygon_centroid_periodic(vpos,verticesOnPolygon,xp::Number,yp::Number)
    output = similar(vpos,length(verticesOnPolygon))
    return compute_polygon_centroid_periodic!(output,vpos,verticesOnPolygon,xp,yp)
end

function compute_polygon_centroid_spherical!(output, vpos, verticesOnPolygon, R::Number)
    @parallel for c in eachindex(verticesOnPolygon)
        @inbounds output[c] = spherical_polygon_centroid(R, vpos,verticesOnPolygon[c])
    end
    return output
end

function compute_polygon_centroid_spherical(vpos, verticesOnPolygon, R::Number)
    output = similar(vpos,length(verticesOnPolygon))
    return compute_polygon_centroid_spherical!(output, vpos, verticesOnPolygon, R)
end

function compute_longitude_periodic(cpos::Vec2DxyArray{T,1}) where {T}
    return zeros(T, length(cpos))
end

function compute_longitude!(longitude::Vector, pos::VecArray)
    px = pos.x
    py = pos.y
    @parallel for c in eachindex(pos)
        @inbounds begin
            longitude[c] = atan(py[c], px[c])
        end
    end
    return longitude
end

function compute_longitude_spherical(cpos::Vec3DArray{T,1}) where {T}
    longitude = Vector{eltype(T)}(undef, length(cpos))
    return compute_longitude!(longitude, cpos)
end

function compute_latitude_periodic(cpos::Vec2DxyArray{T,1}) where {T}
    return zeros(T, length(cpos))
end

function compute_latitude!(latitude::Vector, pos::VecArray)
    px = pos.x
    py = pos.y
    pz = pos.z
    @parallel for c in eachindex(pos)
        @inbounds begin
            x = px[c]
            y = py[c]
            z = pz[c]
            latitude[c] = atan(z / sqrt(x*x + y*y))
        end
    end
    return latitude
end

function compute_latitude_spherical(cpos::Vec3DArray{T,1}) where {T}
    latitude = Vector{eltype(T)}(undef, length(cpos))
    return compute_latitude!(latitude, cpos)
end

function compute_zonalVector_periodic(cpos::Vec2DxyArray{T,1}) where {T}
    return VecArray(x = ones(T, length(cpos)))
end

function compute_zonalVector!(zonalVector::VecArray, pos::VecArray)
    px = pos.x
    py = pos.y
    @parallel for c in eachindex(pos)
        @inbounds begin
            ϕ = atan(py[c], px[c])
            zonalVector[c] = cos(ϕ)*𝐣 - sin(ϕ)𝐢
        end
    end
    return zonalVector
end

function compute_zonalVector_spherical(cpos::Vec3DArray{T,1}) where {T}
    xv = Vector{eltype(T)}(undef, length(cpos))
    yv = Vector{eltype(T)}(undef, length(cpos))
    zonalVector = VecArray(x = xv, y = yv)
    return compute_zonalVector!(zonalVector, cpos)
end

function compute_meridionalVector_periodic(cpos::Vec2DxyArray{T,1}) where {T}
    return VecArray(y = ones(T, length(cpos)))
end

function compute_meridionalVector!(meridionalVector::VecArray, pos::VecArray)
    px = pos.x
    py = pos.y
    pz = pos.z
    @parallel for c in eachindex(pos)
        @inbounds begin
            x = px[c]
            y = py[c]
            z = pz[c]
            ϕ = atan(y, x)
            θ = atan(sqrt(x*x + y*y), z)
            sinϕ = sin(ϕ)
            cosϕ = cos(ϕ)
            sinθ = sin(θ)
            cosθ = cos(θ)
            meridionalVector[c] = sinθ*𝐤 - (cosϕ*cosθ)*𝐢 - (sinϕ*cosθ)*𝐣
        end
    end
    return meridionalVector
end

function compute_meridionalVector_spherical(cpos::Vec3DArray{T,1}) where {T}
    meridionalVector = similar(cpos)
    return compute_meridionalVector!(meridionalVector, cpos)
end

