#########################################################################
# VoronoiMeshes — Periodic Mesh Series (Power of 2)
#
# This example demonstrates how to:
#   - Construct a series of centroidal Voronoi meshes on periodic domains
#   - Vary the number of cells in powers of 2, starting with 20 cells
#   - Export the Voronoi (dual) and Delaunay (triangulation) grids to VTU files
#   - Use the high-level `save` function for automatic format selection
#
# Requirements:
#   - VoronoiMeshes.jl (with extensions)
#   - WriteVTK (for VTU I/O)
#   - DelaunayTriangulation, TensorsLite (optional)
#########################################################################

# Basic modules for mesh construction and manipulation
using VoronoiMeshes
using DelaunayTriangulation
using TensorsLite
using LinearAlgebra

# Required for VTU export (loads extensions automatically)
using WriteVTK

# Create meshes with cell counts in powers of 2, starting with 20
# Powers: 20*2^0=20, 20*2^1=40, 20*2^2=80, 20*2^3=160, etc.
base_cells = 20
num_scales = 12  # Number of different scales to generate

for i in 0:(num_scales-1)
    num_cells = base_cells * (2^i)

    println("Creating periodic mesh with $num_cells cells...")

    # Create a centroidal Voronoi mesh with num_cells cells on a 1×1 periodic domain
    mesh = VoronoiMesh(num_cells, 1.0, 1.0, max_iter=100)

    # Create filenames with the power index
    voronoi_filename = "output/mesh_periodic_p$(i)_voronoi.vtu"
    triangulation_filename = "output/mesh_periodic_p$(i)_triangulation.vtu"
    base_filename = "output/mesh_periodic_p$(i).vtu"

    # --- Export Voronoi (dual) grid to VTU ---
    save_voronoi_to_vtu(voronoi_filename, mesh)
    println("  Saved: $voronoi_filename")

    # --- Export Delaunay triangulation (primal) grid to VTU ---
    save_triangulation_to_vtu(triangulation_filename, mesh)
    println("  Saved: $triangulation_filename")

    # --- Use high-level save function ---
    # This will create two vtu files with _vor and _tri suffixes
    VoronoiMeshes.save(base_filename, mesh)
    println("  Saved: output/mesh_periodic_p$(i)_vor.vtu and output/mesh_periodic_p$(i)_tri.vtu")
    println()
end

println("All meshes created successfully!")
