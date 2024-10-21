function compute_cell_area_periodic!(output,vpos,verticesOnCell,xp::Number,yp::Number)
    @parallel for c in eachindex(verticesOnCell)
        @inbounds output[c] = area(vpos,verticesOnCell[c],xp,yp)
    end
    return output
end

function compute_cell_area_periodic(vpos,verticesOnCell,xp::Number,yp::Number)
    output = Vector{nonzero_eltype(eltype(vpos))}(undef,length(verticesOnCell))
    return compute_cell_area_periodic!(output,vpos,verticesOnCell,xp,yp)
end

compute_cell_area!(output,cells::Cells{false}) = compute_cell_area_periodic!(output, cells.info.diagram.vertices, cells.vertices, cells.x_period, cells.y_period)

compute_cell_area(cells::Cells{false}) = compute_cell_area_periodic(cells.info.diagram.vertices, cells.vertices,cells.x_period, cells.y_period)

function compute_cell_area_spherical!(output, vpos, verticesOnCell, R::Number)
    @parallel for c in eachindex(verticesOnCell)
        @inbounds output[c] = spherical_polygon_area(R, vpos,verticesOnCell[c])
    end
    return output
end

function compute_cell_area_spherical(vpos, verticesOnCell, R::Number)
    output = Vector{nonzero_eltype(eltype(vpos))}(undef, length(verticesOnCell))
    return compute_cell_area_spherical!(output, vpos, verticesOnCell, R)
end

compute_cell_area!(output,cells::Cells{true}) = compute_cell_area_spherical!(output, cells.info.diagram.vertices, cells.vertices, cells.sphere_radius)

compute_cell_area(cells::Cells{true}) = compute_cell_area_spherical(cells.info.diagram.vertices, cells.vertices, cells.sphere_radius)

