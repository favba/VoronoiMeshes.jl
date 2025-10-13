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
mesh = VoronoiMesh(20, 1.0, 1.0)

points = [mesh.cells.position.x'; mesh.cells.position.y'; mesh.cells.position.z']
n = length(mesh.vertices.cells)
triangles = Vector{MeshCell}(undef, n)
Threads.@threads for i in 1:n
    # allocate per-iteration data to avoid races
    idx = collect(mesh.vertices.cells[i])
    triangles[i] = MeshCell(VTKCellTypes.VTK_TRIANGLE, idx)
end
vtk = vtk_grid("test.vtu", points, triangles)
saved_files = close(vtk)

# Print a summary of the mesh
println(mesh)


