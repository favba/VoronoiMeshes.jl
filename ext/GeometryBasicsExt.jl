module GeometryBasicsExt

using VoronoiMeshes
using TensorsLite: Vec, norm
using TensorsLiteGeometry, ImmutableVectors
using GeometryBasics: GeometryBasics, Polygon, Point2f, LineString, Line, TupleView

@static if pkgversion(GeometryBasics) < v"0.5"
    const PolType = Polygon{2, Float32, Point2f, LineString{2, Float32, Point2f, Base.ReinterpretArray{Line{2, Float32}, 1, Tuple{Point2f, Point2f}, TupleView{Tuple{Point2f, Point2f}, 2, 1, Vector{Point2f}}, false}}, Vector{LineString{2, Float32, Point2f, Base.ReinterpretArray{Line{2, Float32}, 1, Tuple{Point2f, Point2f}, TupleView{Tuple{Point2f, Point2f}, 2, 1, Vector{Point2f}}, false}}}}
else
    const PolType = Polygon{2, Float32}
end

function create_cell_polygons_periodic(vert_pos, cell_pos, verticesOnCell, x_period, y_period)

    cell_polygons = Vector{PolType}(undef, length(cell_pos))
    NE = max_length(eltype(verticesOnCell))
    @parallel for i in eachindex(cell_pos)
        @inbounds begin
            cpos = cell_pos[i]
            local_vertices = ImmutableVector{NE, Point2f}()
            for i_v in verticesOnCell[i]
                vpos = closest(cpos, vert_pos[i_v], x_period, y_period)
                local_vertices = @inbounds push(local_vertices, Point2f(vpos.x, vpos.y))
            end
            cell_polygons[i] = Polygon(Array(local_vertices))
        end
    end

    return cell_polygons
end

function VoronoiMeshes.create_cell_polygons(mesh::AbstractVoronoiMesh{false})
    return create_cell_polygons_periodic(mesh.vertices.position, mesh.cells.position, mesh.cells.vertices, mesh.x_period, mesh.y_period)
end

function create_cell_polygons_sphere(vert_lon::Vector{T}, vert_lat, cell_lon, verticesOnCell) where {T<:Number}
    lon_min = minimum(cell_lon)

    if lon_min < zero(T)
        lon_factor = zero(T)
    else
        lon_factor = T(-180)
    end

    x_period = T(360)

    cell_polygons = Vector{PolType}(undef, length(cell_lon))
    NE = max_length(eltype(verticesOnCell))

    @parallel for i in eachindex(cell_lon)
        @inbounds begin
            c_lon = rad2deg(cell_lon[i]) + lon_factor
            local_vertices = ImmutableVector{NE, Point2f}()
            for i_v in verticesOnCell[i]
                vlon_aux = rad2deg(vert_lon[i_v]) + lon_factor
                vlat = rad2deg(vert_lat[i_v])
                vlons = (vlon_aux - x_period, vlon_aux, vlon_aux + x_period)
                _, j = findmin(abs, vlons .- c_lon)
                vlon = vlons[j]
                
                local_vertices = @inbounds push(local_vertices, Point2f(vlon, vlat))
            end
            cell_polygons[i] = Polygon(Array(local_vertices))
        end
    end

    return cell_polygons
end

function VoronoiMeshes.create_cell_polygons(mesh::AbstractVoronoiMesh{true})
    return create_cell_polygons_sphere(mesh.vertices.longitude, mesh.vertices.latitude, mesh.cells.longitude, mesh.cells.vertices)
end

function create_dual_triangles_periodic(vert_pos, cell_pos, cellsOnVertex, x_period, y_period)

    vert_triangles = Vector{PolType}(undef, length(vert_pos))

    @parallel for i in eachindex(vert_pos)
        @inbounds begin
            vpos = vert_pos[i]
            ic1, ic2, ic3 = cellsOnVertex[i]

            cpos1 = closest(vpos, cell_pos[ic1], x_period, y_period)

            cpos2 = closest(vpos, cell_pos[ic2], x_period, y_period)

            cpos3 = closest(vpos, cell_pos[ic3], x_period, y_period)

            vert_triangles[i] = Polygon([Point2f(cpos1.x, cpos1.y), Point2f(cpos2.x, cpos2.y), Point2f(cpos3.x, cpos3.y)])
        end
    end

    return vert_triangles
end

function VoronoiMeshes.create_dual_triangles(mesh::AbstractVoronoiMesh{false})
    return create_dual_triangles_periodic(mesh.vertices.position, mesh.cells.position, mesh.vertices.cells, mesh.x_period, mesh.y_period)
end

function create_edge_quadrilaterals_periodic(edge_pos, vert_pos, cell_pos, verticesOnEdge, cellsOnEdge, x_period, y_period)

    edge_quadrilaterals = Vector{PolType}(undef, length(edge_pos))

    @parallel for i in eachindex(edge_pos)
        @inbounds begin
            epos = edge_pos[i]
            iv1, iv2 = verticesOnEdge[i]
            ic1, ic2 = cellsOnEdge[i]

            cpos1 = closest(epos, cell_pos[ic1], x_period, y_period)

            vpos1 = closest(epos, vert_pos[iv1], x_period, y_period)

            cpos2 = closest(epos, cell_pos[ic2], x_period, y_period)

            vpos2 = closest(epos, vert_pos[iv2], x_period, y_period)

            edge_quadrilaterals[i] = Polygon([Point2f(cpos1.x, cpos1.y), Point2f(vpos1.x, vpos1.y), Point2f(cpos2.x, cpos2.y), Point2f(vpos2.x, vpos2.y)])
        end
    end

    return edge_quadrilaterals
end

function VoronoiMeshes.create_edge_quadrilaterals(mesh::AbstractVoronoiMesh{false})
    return create_edge_quadrilaterals_periodic(mesh.edges.position, mesh.vertices.position, mesh.cells.position, mesh.edges.vertices, mesh.edges.cells, mesh.x_period, mesh.y_period)
end

function create_cell_linesegments_periodic(vert_pos, edge_pos, cell_pos, edgesOnCell, verticesOnEdge, x_period, y_period)
    x = eltype(edge_pos.x)[]
    y = eltype(edge_pos.y)[]
    nEdges = length(verticesOnEdge)
    sizehint!(x, 2 * nEdges)
    sizehint!(y, 2 * nEdges)
    x_periodic = eltype(edge_pos.x)[]
    y_periodic = eltype(edge_pos.y)[]
    sizehint!(x_periodic, nEdges ÷ 2)
    sizehint!(y_periodic, nEdges ÷ 2)

    touched_edges_pos = Set{eltype(edge_pos)}()

    @inbounds for i in eachindex(edgesOnCell)
        c_pos = cell_pos[i]
        edges_ind = edgesOnCell[i]

        for j in edges_ind
            e_pos = edge_pos[j]
            closest_e_pos = closest(c_pos, e_pos, x_period, y_period)
            if !(closest_e_pos in touched_edges_pos)
                p = verticesOnEdge[j]

                is_interior = closest_e_pos == e_pos
                for k in p
                    v_pos = vert_pos[k]
                    closest_v_pos = closest(closest_e_pos, v_pos, x_period, y_period)
                    if (is_interior)
                        push!(x, closest_v_pos.x)
                        push!(y, closest_v_pos.y)
                    else
                        push!(x_periodic, closest_v_pos.x)
                        push!(y_periodic, closest_v_pos.y)
                    end
                end

                push!(touched_edges_pos, closest_e_pos)
            end

        end

    end

    return ((x, y), (x_periodic, y_periodic))
end

function VoronoiMeshes.create_cell_linesegments(mesh::AbstractVoronoiMesh{false})
    return create_cell_linesegments_periodic(mesh.vertices.position, mesh.edges.position, mesh.cells.position, mesh.cells.edges, mesh.edges.vertices, mesh.x_period, mesh.y_period)
end

function create_diagram_linesegments_periodic(vert_pos, cell_pos, verticesOnCell, cellsOnVertex, x_period, y_period)
    x = eltype(vert_pos.x)[]
    y = eltype(vert_pos.y)[]
    nEdges = length(verticesOnCell) + length(cellsOnVertex)
    sizehint!(x, 2 * nEdges)
    sizehint!(y, 2 * nEdges)

    touched_edges_pos = Set{eltype(vert_pos)}()

    @inbounds for i in eachindex(verticesOnCell)
        c_pos = cell_pos[i]
        vertices_ind = verticesOnCell[i]

        L = length(vertices_ind)
        for j in Base.OneTo(L)
            j
            jp1 = j + 1
            j2 = ifelse(jp1 > L, jp1 - L, jp1)
            v1 = vertices_ind[j]
            v2 = vertices_ind[j2]
            v1_pos = closest(c_pos, vert_pos[v1], x_period, y_period)
            v2_pos = closest(c_pos, vert_pos[v2], x_period, y_period)
            e_pos = (v1_pos + v2_pos) / 2
            if !(e_pos in touched_edges_pos)
                push!(x, v1_pos.x)
                push!(x, v2_pos.x)
                push!(y, v1_pos.y)
                push!(y, v2_pos.y)
                push!(touched_edges_pos, e_pos)
            end
        end
    end
    return (x, y)
end

function VoronoiMeshes.create_diagram_linesegments(diagram::VoronoiDiagram{false})
    return create_diagram_linesegments_periodic(diagram.vertices, diagram.generators, diagram.verticesOnCell, diagram.cellsOnVertex, diagram.x_period, diagram.y_period)
end

function create_dual_triangles_linesegments_periodic(vert_pos, edge_pos, cell_pos, edgesOnVertex, cellsOnEdge, x_period, y_period)
    x = eltype(edge_pos.x)[]
    y = eltype(edge_pos.y)[]
    nEdges = length(cellsOnEdge)
    sizehint!(x, 2 * nEdges)
    sizehint!(y, 2 * nEdges)
    x_periodic = eltype(edge_pos.x)[]
    y_periodic = eltype(edge_pos.y)[]
    sizehint!(x_periodic, nEdges ÷ 2)
    sizehint!(y_periodic, nEdges ÷ 2)

    touched_edges_pos = Set{eltype(edge_pos)}()

    @inbounds for i in eachindex(edgesOnVertex)
        v_pos = vert_pos[i]
        edges_ind = edgesOnVertex[i]

        for j in edges_ind
            e_pos = edge_pos[j]
            closest_e_pos = closest(v_pos, e_pos, x_period, y_period)
            if !(closest_e_pos in touched_edges_pos)
                p = cellsOnEdge[j]

                is_interior = closest_e_pos == e_pos
                for k in p
                    c_pos = cell_pos[k]
                    closest_c_pos = closest(closest_e_pos, c_pos, x_period, y_period)
                    if (is_interior)
                        push!(x, closest_c_pos.x)
                        push!(y, closest_c_pos.y)
                    else
                        push!(x_periodic, closest_c_pos.x)
                        push!(y_periodic, closest_c_pos.y)
                    end
                end

                push!(touched_edges_pos, closest_e_pos)
            end

        end

    end

    return ((x, y), (x_periodic, y_periodic))
end

function VoronoiMeshes.create_dual_triangles_linesegments(mesh::AbstractVoronoiMesh{false})
    return create_dual_triangles_linesegments_periodic(mesh.vertices.position, mesh.edges.position, mesh.cells.position, mesh.vertices.edges, mesh.edges.cells, mesh.x_period, mesh.y_period)
end

end

