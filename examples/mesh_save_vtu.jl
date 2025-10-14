#########################################################################
# VoronoiMeshes — vtu save
#
#########################################################################
using VoronoiMeshes
using DelaunayTriangulation
using TensorsLite
using LinearAlgebra

using GLMakie  # Optional, for plotting if Makie is installed
using NCDatasets # Optional, for saving/loading meshes in NetCDF format
using ReadVTK  # For loading meshes in VTU format
using WriteVTK # For saving meshes in VTU format

# Create a centroidal Voronoi mesh with 20 cells on a 1×1 periodic domain
mesh = VoronoiMesh(40, 1.0, 1.0)

points = [mesh.cells.position.x'; mesh.cells.position.y'; mesh.cells.position.z']
nverts = length(mesh.vertices.cells)
triangles = Vector{MeshCell}(undef, nverts)
for i in 1:nverts
    # allocate per-iteration data to avoid races
    idx = collect(mesh.vertices.cells[i])
    triangles[i] = MeshCell(VTKCellTypes.VTK_TRIANGLE, idx)
end
vtk = vtk_grid("mesh_triangular.vtu", points, triangles)
saved_files = close(vtk)


# Print a summary of the mesh
println(mesh)

# --- write dual (Voronoi) grid to VTU ---

points_vertices = [mesh.vertices.position.x'; mesh.vertices.position.y'; mesh.vertices.position.z']

ncells = mesh.cells.n
voronoi_cells = Vector{MeshCell}(undef, ncells)
for i in 1:ncells
    idx = collect(mesh.cells.vertices[i])
    # Voronoi cells are polygons; use the VTK polygon cell type instead of
    # constructing internal PolyData objects (those are opaque and have no fields).
    voronoi_cells[i] = MeshCell(VTKCellTypes.VTK_POLYGON, idx)
end

vtk2 = vtk_grid("mesh_voronoi.vtu", points_vertices, voronoi_cells)
saved_files2 = close(vtk2)


# --- write triangular grid on a edge by edge way ---

#points_vertices = [mesh.vertices.position.x'; mesh.vertices.position.y'; mesh.vertices.position.z']


let
    vor_poly = create_cell_polygons(mesh)
    n_polys = length(vor_poly)
    points_mat = Matrix{Float64}(undef, 2, n_polys*6)  # preallocate assuming average 6 vertices per polygon
    voronoi_cells = Vector{MeshCell}(undef, n_polys)
    ivtx = Int64(1)
    for i in 1:n_polys
        ivtx_start = ivtx
        for pt in vor_poly[i]
            points_mat[1, ivtx] = pt[1]
            points_mat[2, ivtx] = pt[2]
            ivtx += 1
        end
        links = collect(ivtx_start:ivtx - 1)
        voronoi_cells[i] = MeshCell(VTKCellTypes.VTK_POLYGON, links)
    end
    println(points_mat)
    vtk2 = vtk_grid("mesh_periodic.vtu", points_mat, voronoi_cells)
    saved_files2 = close(vtk2)
end

