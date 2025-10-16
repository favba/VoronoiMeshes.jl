
#########################################################################
# VoronoiMeshes — VTU save/read example
#
# This example demonstrates how to:
#   - Construct a centroidal Voronoi mesh on a periodic domain
#   - Export the Voronoi (dual) and Delaunay (triangulation) grids to VTU files
#   - Use the high-level `save` function for automatic format selection
#
# Requirements:
#   - VoronoiMeshes.jl (with extensions)
#   - WriteVTK, VTKBase, ReadVTK (for VTU I/O)
#   - DelaunayTriangulation, TensorsLite, GLMakie, NCDatasets (optional)
#########################################################################

# Basic modules for mesh construction and manipulation
using VoronoiMeshes
using DelaunayTriangulation
using TensorsLite
using LinearAlgebra

# Optional: for plotting and NetCDF I/O
using GLMakie      # For plotting (if Makie is installed)
using NCDatasets   # For saving/loading meshes in NetCDF format

# Required for VTU export/import (loads extensions automatically)
using WriteVTK     # For saving meshes in VTU format

# Create a centroidal Voronoi mesh with 40 cells on a 1×1 periodic domain
mesh = VoronoiMesh(40, 1.0, 1.0)

# --- Export Voronoi (dual) grid to VTU ---
# This will handle periodic ghost vertices automatically
save_voronoi_to_vtu("mesh_voronoi.vtu", mesh)

# --- Export Delaunay triangulation (primal) grid to VTU ---
save_triangulation_to_vtu("mesh_triangulation.vtu", mesh)

# --- Use high-level save function (auto-selects format by extension) ---
# This will call the appropriate extension based on the file extension
VoronoiMeshes.save("mesh.vtu", mesh)

mesh =  VoronoiMesh("mesh.vtu")
