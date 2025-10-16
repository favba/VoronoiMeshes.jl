

# ReadVTKExt: Import a VoronoiMesh meshes from VTU (VTK Unstructured Grid) format.
# - Provides functions to import Voronoi and Delaunay meshes with periodic ghost vertices.
# - Uses ReadVTK for file input
module ReadVTKExt

using VoronoiMeshes, TensorsLite, TensorsLiteGeometry, Zeros, SmallCollections, LinearAlgebra
import VoronoiMeshes: save_triangulation_to_vtu, save_voronoi_to_vtu

using PrecompileTools

using ReadVTK # For saving meshes in VTU format
using ReadVTK.VTKBase


#Wrapper to load a VoronoiMesh from VTU files with a single filename
function VoronoiMeshes.VoronoiMesh(filename::String)
    # If a single filename was given, search for the 2 associated grids (vor and tri)
    println("Attempt to read VTU file with base name: ", filename)
    name, ext = Base.Filesystem.splitext(filename)
    if ext == ".vtu" #Save to VTU using VTKExt
        name_vor = name*"_vor"*ext
        name_tri = name*"_tri"*ext
    else
        error("Unsupported file extension: $filename")
    end
    return VoronoiMeshes.VoronoiMesh(name_vor, name_tri)
end

# Function to load a VoronoiMesh from VTU files (2 files must be given)
function VoronoiMeshes.VoronoiMesh(filename_vor::String, filename_tri::String)
    if isfile(filename_vor) && isfile(filename_tri)
            println("Loading files: ", filename_vor, " and ", filename_tri)
        else
            error("Couldn't find files: ", filename_vor, " and ", filename_tri)
        end
    vtk_tri = VTKFile(filename_tri)
    vtk_vor = VTKFile(filename_vor)
    println("VTU triangulation file read: ", vtk_tri)

    # this is what we need to reconstruct the VoronoiMesh
    # generators::Vec2DxyArray{TF, 1}
    # vertices::Vec2DxyArray{TF, 1}
    # verticesOnCell::SmVecArray{maxEdges, TI, 1}
    # cellsOnVertex::Vector{FixedVector{3, TI}}
    # meshDensity::Vector{TF}
    # x_period::TF
    # y_period::TF
    cell_data_tri = get_cell_data(vtk_tri)
    point_data_tri = get_point_data(vtk_tri)
    vtk_cells_tri = get_cells(vtk_tri)
    vtk_points_tri = get_points(vtk_tri)

    return 0
end

end # module VTKExt
