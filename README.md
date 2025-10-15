# VoronoiMeshes

VoronoiMeshes.jl provides tools to create, inspect and save biperiodic planar and spherical Voronoi meshes (including centroidal/Lloyd methods and utilities to work with their dual Delaunay triangulations).

Main developer: Felipe A. V. de Bragança Alves <favbalves@gmail.com>

![Bi-periodic planar mesh creation](./assets/mesh_creation_15fps_748x551.avif)

## Install Guide 

Assuming Julia already installed (for instance via juliaup). 

The package, and the dependencies, will be ported into registered packages soon. 

```julia
import Pkg
Pkg.activate(temp=true)
Pkg.add("DelaunayTriangulation")
Pkg.add(url="https://github.com/favba/TensorsLite.jl.git")
Pkg.add(url="https://github.com/favba/TensorsLiteGeometry.jl.git")
Pkg.add(url="https://github.com/favba/VoronoiMeshes.jl.git")
```


## Developer's Guide 

Assuming Julia already installed (for instance via juliaup). 

The package, and the dependencies, will be ported into registered packages soon. For current developers, the enviroment can be constructed in the following way.

Assuming a comom folder for the packages and the dependencies:
```bash
git clone git@github.com:favba/VoronoiMeshes.jl.git
git clone git@github.com:favba/TensorsLite.jl.git
git clone git@github.com:favba/TensorsLiteGeometry.jl.git
```

Create project (manisfests) for VoronoiMeshes

```
cd TensorsLite.jl
julia --project=. -e 'import Pkg; Pkg.instantiate(); Pkg.status()'
cd ../TensorsLiteGeometry.jl
julia --project=. -e 'import Pkg; Pkg.develop(path="../TensorsLite.jl"); Pkg.instantiate(); Pkg.status()'
cd ../VoronoiMeshes.jl
julia --project=. -e 'import Pkg; Pkg.develop([Pkg.PackageSpec(path="../TensorsLite.jl"), Pkg.PackageSpec(path="../TensorsLiteGeometry.jl")]); Pkg.instantiate(); Pkg.status()'
```

Create project (manifests) for 'test'
```
cd test

julia --project=. -e 'import Pkg; Pkg.develop([Pkg.PackageSpec(path="../../TensorsLite.jl"), Pkg.PackageSpec(path="../../TensorsLiteGeometry.jl"), Pkg.PackageSpec(path="../../VoronoiMeshes.jl")]); Pkg.add("DelaunayTriangulation"); Pkg.add("GeometryBasics"); Pkg.add("GLMakie"); Pkg.add("LinearAlgebra"); Pkg.add("VTKBase"); Pkg.add("ReadVTK"); Pkg.add("WriteVTK"); Pkg.add("NCDatasets"); Pkg.add("SmallCollections"); Pkg.add("Test"); Pkg.instantiate(); Pkg.status()'
```


## Quick usage guide


Creating a simple centroidal biperiodic planar Voronoi mesh from N random generators:

```julia 
import Pkg
Pkg.add("DelaunayTriangulation") 
using VoronoiMeshes, DelaunayTriangulation

# Create a centroidal Voronoi mesh with 200 cells on a 1×1 periodic domain
mesh = VoronoiMesh(20, 1.0, 1.0)
println(mesh)

#Mesh plotting
using GLMakie
plotmesh(mesh)
plotdualmesh!(mesh) # Plot on top of the previous plot
```

Create a mesh from a given set of generator points (VecArray of x,y coordinates):

```julia
using VoronoiMeshes, DelaunayTriangulation, TensorsLite

generators = VecArray(x = rand(10), y = rand(10))
mesh = VoronoiMesh(generators, 1.0, 1.0)
println(mesh)
```

Create a regular hexagonal planar mesh (utility provided in the Delaunay extension):

```julia
using DelaunayTriangulation, VoronoiMeshes

# dx ~ target cell spacing
hexmesh = create_planar_hex_mesh(1.0, 1.0, 0.05)
```

