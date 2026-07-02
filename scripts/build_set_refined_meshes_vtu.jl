# build_set_refined_meshes_vtu.jl
#
# Builds a series of centroidal Voronoi meshes with increasing cell counts,
# each generated independently from random initial points and converged via
# Lloyd's algorithm. Cell counts grow as base_cells * 2^i.
#
# Usage: run interactively or adjust base_cells / num_scales / ini_scale below.
#
# Requirements: VoronoiMeshes.jl, WriteVTK, DelaunayTriangulation, GLMakie

# Basic modules for mesh construction and manipulation
using VoronoiMeshes
using DelaunayTriangulation
using TensorsLite
using LinearAlgebra

# Required for VTU export (loads extensions automatically)
using WriteVTK
using GLMakie  # swap for CairoMakie on headless servers

# Create meshes with cell counts in powers of 2, starting with 20
# Powers: 20*2^0=20, 20*2^1=40, 20*2^2=80, 20*2^3=160, etc.
base_cells = 20
num_scales = 13  # Number of different scales to generate
ini_scale = 0     # Starting power index (0 corresponds to base_cells)

function save_mesh_png(filename, mesh, label)
    fig = Figure(size=(700, 700))
    ax = Axis(fig[1, 1], title=label, aspect=DataAspect())
    plotdualmesh!(ax, mesh)
    plotmesh!(ax, mesh)
    hidedecorations!(ax)
    GLMakie.save(filename, fig)
    return nothing
end

for i in ini_scale:(num_scales-1)
    num_cells = base_cells * (2^i)

    println("Creating periodic mesh with $num_cells cells...")

    # Create a centroidal Voronoi mesh with num_cells cells on a 1×1 periodic domain
    mesh = VoronoiMesh(num_cells, 1.0, 1.0, rtol=1e-4, max_iter=10000) #, max_iter=10)

    # Create filenames with the power index
    voronoi_filename = "output/mesh_periodic_refined_p$(i)_voronoi.vtu"
    triangulation_filename = "output/mesh_periodic_refined_p$(i)_triangulation.vtu"
    base_filename = "output/mesh_periodic_refined_p$(i).vtu"

    # --- Export Voronoi (dual) grid to VTU ---
    save_voronoi_to_vtu(voronoi_filename, mesh)
    println("  Saved: $voronoi_filename")

    # --- Export Delaunay triangulation (primal) grid to VTU ---
    save_triangulation_to_vtu(triangulation_filename, mesh)
    println("  Saved: $triangulation_filename")

    # --- Use high-level save function ---
    # This will create two vtu files with _vor and _tri suffixes
    VoronoiMeshes.save(base_filename, mesh)
    println("  Saved: output/mesh_periodic_refined_p$(i)_vor.vtu and output/mesh_periodic_refined_p$(i)_tri.vtu")

    # --- PNG plot (Voronoi + triangulation overlaid) ---
    label = "p$i — $(num_cells) cells"
    save_mesh_png("output/mesh_periodic_refined_p$(i).png", mesh, label)
    println("  Saved: output/mesh_periodic_refined_p$(i).png")
    println()
end

println("All meshes created successfully!")
