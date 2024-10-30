compute_vertex_area!(output, vertices::Vertices{false}) = compute_polygon_area_periodic!(output, vertices.info.diagram.generators, vertices.cells, vertices.x_period, vertices.y_period)
compute_vertex_area!(output, vertices::Vertices{true}) = compute_polygon_area_spherical!(output, vertices.info.diagram.generators, vertices.cells, vertices.sphere_radius)
compute_vertex_area(vertices::Vertices) = compute_vertex_area!(similar(vertices.position.x), vertices)

compute_vertex_centroid!(output, vertices::Vertices{false}) = compute_polygon_centroid_periodic!(output, vertices.info.diagram.generators, vertices.cells, vertices.x_period, vertices.y_period)
compute_vertex_centroid!(output, vertices::Vertices{true}) = compute_polygon_centroid_spherical!(output, vertices.info.diagram.generators, vertices.cells, vertices.sphere_radius)
compute_vertex_centroid(vertices::Vertices) = compute_vertex_centroid!(similar(vertices.position), vertices)

compute_vertex_longitude!(output, vertices::Vertices{false}) = fill!(output, zero(float_type(vertices)))
compute_vertex_longitude!(output, vertices::Vertices{true}) = compute_longitude!(output, vertices.position)
compute_vertex_longitude(vertices::Vertices) = compute_vertex_longitude!(similar(vertices.position.x), vertices)

compute_vertex_latitude!(output, vertices::Vertices{false}) = fill!(output, zero(float_type(vertices)))
compute_vertex_latitude!(output, vertices::Vertices{true}) = compute_latitude!(output, vertices.position)
compute_vertex_latitude(vertices::Vertices) = compute_vertex_latitude!(similar(vertices.position.x), vertices)

function compute_kite_areas_periodic!(output, cpos, vpos, cellsOnVertex, xp, yp)

    @parallel for v in eachindex(cellsOnVertex)
        @inbounds begin
            v_pos = vpos[v]
            c1, c2, c3 = cellsOnVertex[v]
            c1_pos = closest(v_pos, cpos[c1], xp, yp)
            c2_pos = closest(v_pos, cpos[c2], xp, yp)
            c3_pos = closest(v_pos, cpos[c3], xp, yp)

            #Those should be the edges positions
            c12_pos = 0.5 * (c1_pos + c2_pos)
            c23_pos = 0.5 * (c2_pos + c3_pos)
            c31_pos = 0.5 * (c3_pos + c1_pos)

            a1 = area(c1_pos, c12_pos, v_pos, c31_pos)
            a2 = area(c12_pos, c2_pos, c23_pos, v_pos)
            a3 = area(c23_pos, c3_pos, c31_pos, v_pos)
            output[v] = (a1, a2, a3)
        end
    end
    return output
end

function compute_kite_areas_spherical!(output, cpos, vpos, cellsOnVertex, R)

    @parallel for v in eachindex(cellsOnVertex)
        @inbounds begin
            v_pos = vpos[v]
            c1, c2, c3 = cellsOnVertex[v]
            c1_pos = cpos[c1]
            c2_pos = cpos[c2]
            c3_pos = cpos[c3]

            #Those should be the edges positions
            c12_pos = arc_midpoint(R, c1_pos, c2_pos)
            c23_pos = arc_midpoint(R, c2_pos, c3_pos)
            c31_pos = arc_midpoint(R, c3_pos, c1_pos)

            a1 = spherical_polygon_area(R, c1_pos, c12_pos, v_pos, c31_pos)
            a2 = spherical_polygon_area(R, c12_pos, c2_pos, c23_pos, v_pos)
            a3 = spherical_polygon_area(R, c23_pos, c3_pos, c31_pos, v_pos)
            output[v] = (a1, a2, a3)
        end
    end
    return output
end

compute_vertex_kiteAreas!(output, vertices::Vertices{false}) = compute_kite_areas_periodic!(
    output, vertices.info.diagram.generators, vertices.position,
    vertices.cells, vertices.x_period, vertices.y_period
)
compute_vertex_kiteAreas!(output, vertices::Vertices{true}) = compute_kite_areas_spherical!(
    output, vertices.info.diagram.generators, vertices.position,
    vertices.cells, vertices.sphere_radius
)
compute_vertex_kiteAreas(vertices::Vertices) = compute_vertex_kiteAreas!(Vector{NTuple{3, float_type(vertices)}}(undef, vertices.n), vertices)
