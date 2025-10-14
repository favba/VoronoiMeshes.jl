module WriteVTKExt

using VoronoiMeshes, TensorsLite, Zeros, SmallCollections
import VoronoiMeshes: save_to_vtk, save_to_vtk!, copy_matrix_to_fixedvector_vector!
using PrecompileTools
using WriteVTK
using VTKBase  # if required for VTKCellTypes
using Threads
using LinearAlgebra

points = [mesh.cells.position.x'; mesh.cells.position.y'; mesh.cells.position.z']
n = length(mesh.vertices.cells)
triangles = Vector{MeshCell}(undef, n)
Threads.@threads for i in 1:n
           # allocate per-iteration data to avoid races
           idx = collect(mesh.vertices.cells[i])
           triangles[i] = MeshCell(VTKCellTypes.VTK_TRIANGLE, idx)
       end

vtk = vtk_grid("mesh_primal.vtu", points, triangles)
saved_files = close(vtk)

end # module WriteVTKExt
