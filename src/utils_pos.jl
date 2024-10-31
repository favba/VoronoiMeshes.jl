function graph_partition(cellsOnCell::AbstractVector{<:ImmutableVector},nEdges::Integer)
    nCells = length(cellsOnCell)
    outIO = IOBuffer()
    println(outIO, nCells,' ',nEdges)

    @inbounds for i in eachindex(cellsOnCell)
        map(x->print(outIO,x,' '),cellsOnCell[i])
        skip(outIO,-1)
        println(outIO)
    end
    seekstart(outIO)
    return outIO
end

graph_partition(cells::Cells,edges::Edges) = graph_partition(cells.cells,edges.n)

graph_partition(mesh::VoronoiMesh) = graph_partition(mesh.cells,mesh.edges)

function find_obtuse_triangles_periodic(cpos,cellsOnVertex,xp::Number,yp::Number)
    r = Int[]
    lk = ReentrantLock()
    @parallel for v in eachindex(cellsOnVertex)
        @inbounds begin
        c1,c2,c3 = cellsOnVertex[v]
        c1pos = cpos[c1]
        c2pos = closest(c1pos,cpos[c2],xp,yp)
        c3pos = closest(c1pos,cpos[c3],xp,yp)
        if is_obtuse(c1pos,c2pos,c3pos)
            lock(lk) do
                push!(r,v)
            end
        end
        end #inbounds
    end
    return r
end

find_obtuse_triangles(vertices::Vertices{false},cells::Cells{false},xp::Number,yp::Number) = find_obtuse_triangles_periodic(cells.position,vertices.cells,xp,yp)

find_obtuse_triangles(mesh::VoronoiMesh{false}) = find_obtuse_triangles(mesh.vertices,mesh.cells,mesh.x_period, mesh.y_period)

"""
    periodic_edges_mask(dcEdge::Vector,cellsOnEdge,cell_positions)
    periodic_edges_mask(mesh::VoronoiMesh)

Returns a BitArray that masks periodic edges (interior edges = true, boundary edges = false).
"""
function periodic_edges_mask(dc,cellsOnEdge,c_pos)
    mask = BitArray(undef,length(dc))
    @inbounds for e in eachindex(cellsOnEdge)
        dc_e = dc[e]
        c1,c2 = cellsOnEdge[e]
        dc_2 = norm(c_pos[c2] - c_pos[c1])
        mask[e] = dc_e ≈ dc_2
    end
    return mask
end

periodic_edges_mask(mesh::VoronoiMesh{false}) = periodic_edges_mask(mesh.edges.cellsDistance, mesh.edges.cells, mesh.cells.position)

function periodic_edges_mask(mesh::VoronoiMesh{true})
    a = BitArray(mesh.edges.n)
    a .= true
    return a
end

"""
    periodic_vertices_mask(mesh::VoronoiMesh)

Returns a BitArray that masks vertices that belong to a periodic edge (interior vertices = true, boundary vertices = false).
"""
function periodic_vertices_mask(edgesOnVertex, mask_edges)
    mask = BitArray(undef,length(edgesOnVertex))
    @inbounds for v in eachindex(edgesOnVertex)
        e1, e2, e3 = edgesOnVertex[v]
        mask[v] = mask_edges[e1] && mask_edges[e2] && mask_edges[e3]
    end
    return mask
end

periodic_vertices_mask(mesh::VoronoiMesh{false}) = periodic_vertices_mask(mesh.vertices.edges, periodic_edges_mask(mesh))

function periodic_vertices_mask(mesh::VoronoiMesh{true})
    a = BitArray(mesh.vertices.n)
    a .= true
    return a
end
