module GeometryBasicsExt

using VoronoiMeshes
using TensorsLite: Vec, norm
using TensorsLiteGeometry, ImmutableVectors
using GeometryBasics: GeometryBasics, Polygon, Point2f, LineString, Line, TupleView

const PolTypeIm{N} = Polygon{2,Float32,Point2f,LineString{2,Float32,Point2f,Base.ReinterpretArray{GeometryBasics.Ngon{2,Float32,2,Point2f},1,Tuple{Point2f,Point2f},TupleView{Tuple{Point2f,Point2f},2,1,ImmutableVector{N, Point2f}},false}},Array{LineString{2,Float32,Point2f,Base.ReinterpretArray{GeometryBasics.Ngon{2,Float32,2,Point2f},1,Tuple{Point2f,Point2f},TupleView{Tuple{Point2f,Point2f},2,1,ImmutableVector{N, Point2f}},false}},1}} 

function VoronoiMeshes.create_cells_polygons_periodic(vert_pos,cell_pos,verticesOnCell,x_period,y_period)

    cell_polygons = Vector{PolTypeIm{max_length(eltype(verticesOnCell))}}(undef,length(cell_pos))
    NE = max_length(eltype(verticesOnCell))
    @parallel for i in eachindex(cell_pos)
        @inbounds begin
            cpos = cell_pos[i]
            local_vertices = ImmutableVector{NE, Point2f}()
            for i_v in verticesOnCell[i]
                vpos = closest(cpos,vert_pos[i_v],x_period,y_period)
                local_vertices = @inbounds push(local_vertices, Point2f(vpos.x,vpos.y))
            end
            cell_polygons[i] = Polygon(local_vertices)
        end
    end

    return cell_polygons
end

function VoronoiMeshes.create_cells_polygons(mesh::VoronoiMesh{false})
    return VoronoiMeshes.create_cells_polygons_periodic(mesh.vertices.position, mesh.cells.position, mesh.cells.vertices, mesh.x_period, mesh.y_period) 
end

const PolTypeTriangle = PolTypeIm{3}

function VoronoiMeshes.create_dual_triangles_periodic(vert_pos,cell_pos,cellsOnVertex,x_period,y_period)

    vert_triangles = Vector{PolTypeTriangle}(undef,length(vert_pos))

    @parallel for i in eachindex(vert_pos)
        @inbounds begin
            vpos = vert_pos[i]
            ic1,ic2,ic3 = cellsOnVertex[i]

            cpos1 = closest(vpos,cell_pos[ic1],x_period,y_period)

            cpos2 = closest(vpos,cell_pos[ic2],x_period,y_period)

            cpos3 = closest(vpos,cell_pos[ic3],x_period,y_period)

            vert_triangles[i] = Polygon(ImmutableVector{3}(Point2f(cpos1.x,cpos1.y),Point2f(cpos2.x,cpos2.y),Point2f(cpos3.x,cpos3.y)))
        end
    end

    return vert_triangles
end

function VoronoiMeshes.create_dual_triangles(mesh::VoronoiMesh{false})
    return VoronoiMeshes.create_dual_triangles_periodic(mesh.vertices.position, mesh.cells.position, mesh.vertices.cells, mesh.x_period, mesh.y_period) 
end

function VoronoiMeshes.create_edge_quadrilaterals_periodic(edge_pos,vert_pos,cell_pos,verticesOnEdge,cellsOnEdge,x_period,y_period)

    edge_quadrilaterals = Vector{PolTypeIm{4}}(undef,length(edge_pos))

    @parallel for i in eachindex(edge_pos)
        @inbounds begin
            epos = edge_pos[i]
            iv1,iv2 = verticesOnEdge[i]
            ic1,ic2 = cellsOnEdge[i]

            cpos1 = closest(epos,cell_pos[ic1],x_period,y_period)

            vpos1 = closest(epos,vert_pos[iv1],x_period,y_period)

            cpos2 = closest(epos,cell_pos[ic2],x_period,y_period)

            vpos2 = closest(epos,vert_pos[iv2],x_period,y_period)

            edge_quadrilaterals[i] = Polygon(ImmutableVector{4}(Point2f(cpos1.x,cpos1.y),Point2f(vpos1.x,vpos1.y),Point2f(cpos2.x,cpos2.y),Point2f(vpos2.x,vpos2.y)))
        end
    end

    return edge_quadrilaterals
end

function VoronoiMeshes.create_edge_quadrilaterals(mesh::VoronoiMesh{false})
    return VoronoiMeshes.create_edge_quadrilaterals_periodic(mesh.edges.position, mesh.vertices.position, mesh.cells.position, mesh.edges.vertices, mesh.edges.cells, mesh.x_period, mesh.y_period) 
end

function VoronoiMeshes.create_cell_linesegments_periodic(vert_pos,edge_pos,cell_pos,edgesOnCell,verticesOnEdge,x_period,y_period)
    x = eltype(edge_pos.x)[]
    y = eltype(edge_pos.y)[]
    nEdges = length(verticesOnEdge)
    sizehint!(x,2*nEdges)
    sizehint!(y,2*nEdges)
    x_periodic = eltype(edge_pos.x)[]
    y_periodic = eltype(edge_pos.y)[]
    sizehint!(x_periodic,nEdges÷2) 
    sizehint!(y_periodic,nEdges÷2)

    touched_edges_pos = Set{eltype(edge_pos)}()

    @inbounds for i in eachindex(edgesOnCell)
        c_pos = cell_pos[i]
        edges_ind = edgesOnCell[i]

        for j in edges_ind
            e_pos=edge_pos[j]
            closest_e_pos = closest(c_pos,e_pos,x_period,y_period)
            if !(closest_e_pos in touched_edges_pos)
                p = verticesOnEdge[j]

                for k in p
                    v_pos = vert_pos[k]
                    closest_v_pos = closest(closest_e_pos,v_pos,x_period,y_period)
                    if (closest_e_pos == e_pos)
                        push!(x, closest_v_pos.x)
                        push!(y, closest_v_pos.y)
                    else
                        push!(x_periodic, closest_v_pos.x)
                        push!(y_periodic, closest_v_pos.y)
                    end
                end

                push!(touched_edges_pos,closest_e_pos)
            end

        end

    end

    return ((x, y), (x_periodic, y_periodic))
end

function VoronoiMeshes.create_cell_linesegments(mesh::VoronoiMesh{false})
    return VoronoiMeshes.create_cell_linesegments_periodic(mesh.vertices.position, mesh.edges.position, mesh.cells.position, mesh.cells.edges, mesh.edges.vertices, mesh.x_period, mesh.y_period) 
end

end
