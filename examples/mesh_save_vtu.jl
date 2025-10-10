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

# Print a summary of the mesh
println(mesh)


