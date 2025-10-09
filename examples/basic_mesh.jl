# Example: create a simple centroidal Voronoi mesh
using VoronoiMeshes
using DelaunayTriangulation
using TensorsLite

# Create a centroidal Voronoi mesh with 200 cells on a 1×1 periodic domain
mesh1 = VoronoiMesh(20, 1.0, 1.0)

println(mesh1)

# Create a mesh from a given set of generator points (VecArray of x,y coordinates):

generators = VecArray(x = rand(10), y = rand(10))
mesh2 = VoronoiMesh(generators, 1.0, 1.0)
println(mesh2)

# Create a regular hexagonal planar mesh (utility provided in the Delaunay extension):
# dx ~ target cell spacing

hexmesh = create_planar_hex_mesh(1.0, 1.0, 0.05)
println(hexmesh)