
# VTKExt: Extension for VoronoiMeshes to export meshes to VTU format using WriteVTK/VTKBase
module VTKExt

using VoronoiMeshes, TensorsLite, TensorsLiteGeometry, Zeros, SmallCollections, LinearAlgebra
import VoronoiMeshes: save_triangulation_to_vtu, save_voronoi_to_vtu

using PrecompileTools

using VTKBase  # if required for VTKCellTypes
using ReadVTK  # For loading meshes in VTU format
using WriteVTK # For saving meshes in VTU format


"""
    create_ghost_periodic_voronoi_vertices(mesh::AbstractVoronoiMesh{false})

Create Voronoi cell vertex and topology arrays, adding ghost vertices for periodic boundaries.
Returns (vertices_with_ghosts, verticesOnPolygon_with_ghosts, n_ghosts).
"""
function create_ghost_periodic_voronoi_vertices(mesh::AbstractVoronoiMesh{false})
    # Mesh info
    vert_pos = mesh.vertices.position # Voronoi cell vertices
    polygon_pos = mesh.cells.position # Voronoi cell centers
    verticesOnPolygon = mesh.cells.vertices # Indexes of vertices per Voronoi cell
    return create_ghost_periodic_points(vert_pos, polygon_pos, verticesOnPolygon, mesh.x_period, mesh.y_period)
end



"""
    create_ghost_periodic_triangulation_vertices(mesh::AbstractVoronoiMesh{false})

Create Delaunay triangulation vertex and topology arrays, adding ghost vertices for periodic boundaries.
Returns (vertices_with_ghosts, verticesOnPolygon_with_ghosts, n_ghosts).
"""
function create_ghost_periodic_triangulation_vertices(mesh::AbstractVoronoiMesh{false})
    # Mesh info
    vert_pos = mesh.cells.position # Triangle circumcenters
    polygon_pos = mesh.vertices.position # Triangle vertices
    verticesOnPolygon = mesh.vertices.cells # Indexes of vertices per triangle
    return create_ghost_periodic_points(vert_pos, polygon_pos, verticesOnPolygon, mesh.x_period, mesh.y_period)
end


"""
    create_ghost_periodic_points(vert_pos, polygon_pos, verticesOnPolygon, x_period, y_period)

Helper to add ghost vertices for periodic boundaries to a set of polygons.
Returns (vertices_with_ghosts, verticesOnPolygon_with_ghosts, n_ghosts).
"""
function create_ghost_periodic_points(vert_pos, polygon_pos, verticesOnPolygon, x_period, y_period)
    # Ghost info
    vertices_with_ghosts = copy(vert_pos)
    ivertices_with_ghosts = length(vert_pos) + 1
    verticesOnPolygon_with_ghosts = [Vector(v) for v in verticesOnPolygon]
    for i in eachindex(polygon_pos)
        ppos = polygon_pos[i]
        for j in eachindex(verticesOnPolygon[i])
            vpos = closest(ppos, vert_pos[verticesOnPolygon[i][j]], x_period, y_period)
            if norm(vpos-vert_pos[verticesOnPolygon[i][j]]) > 0.4*min(x_period, y_period)
                push!(vertices_with_ghosts, vpos)
                verticesOnPolygon_with_ghosts[i][j] = ivertices_with_ghosts
                ivertices_with_ghosts += 1
            end
        end
    end
    n_ghosts = length(vertices_with_ghosts) - length(vert_pos)
    println("Added $n_ghosts ghost vertices to handle periodicity in mesh")
    return vertices_with_ghosts, verticesOnPolygon_with_ghosts, n_ghosts
end



"""
    save_voronoi_to_vtu(file_name::String, mesh::AbstractVoronoiMesh{false})

Export the Voronoi mesh to a VTU file, handling periodic ghost vertices.
"""
function save_voronoi_to_vtu(file_name::String, mesh::AbstractVoronoiMesh{false})
    vertices_with_ghosts, verticesOnPolygon_with_ghosts, _ =
        create_ghost_periodic_voronoi_vertices(mesh)
    n_polys = mesh.cells.n
    points_mat = hcat((v[1:2] for v in vertices_with_ghosts)...)        
    cells = Vector{MeshCell}(undef, n_polys)
    for i in 1:n_polys
        links = collect(verticesOnPolygon_with_ghosts[i])
        cells[i] = MeshCell(VTKCellTypes.VTK_POLYGON, links)
    end
    vtk = vtk_grid(file_name, points_mat, cells)
    saved_file = close(vtk)
    return saved_file
end


"""
    save_triangulation_to_vtu(file_name::String, mesh::AbstractVoronoiMesh{false})

Export the Delaunay triangulation (dual mesh) to a VTU file, handling periodic ghost vertices.
"""
function save_triangulation_to_vtu(file_name::String, mesh::AbstractVoronoiMesh{false})
    vertices_with_ghosts, verticesOnPolygon_with_ghosts, _ =
        create_ghost_periodic_triangulation_vertices(mesh)
    n_polys = mesh.vertices.n
    points_mat = hcat((v[1:2] for v in vertices_with_ghosts)...)        
    cells = Vector{MeshCell}(undef, n_polys)
    for i in 1:n_polys
        links = collect(verticesOnPolygon_with_ghosts[i])
        cells[i] = MeshCell(VTKCellTypes.VTK_TRIANGLE, links)
    end
    vtk = vtk_grid(file_name, points_mat, cells)
    saved_file = close(vtk)
    return saved_file
end

end # module VTKExt
