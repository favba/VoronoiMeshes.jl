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
julia --project=. -e 'import Pkg; Pkg.develop([Pkg.PackageSpec(path="../../TensorsLite.jl"), Pkg.PackageSpec(path="../../TensorsLiteGeometry.jl"), Pkg.PackageSpec(path="../../VoronoiMeshes.jl")]); Pkg.add("DelaunayTriangulation"); Pkg.add("GeometryBasics"); Pkg.add("GLMakie"); Pkg.add("LinearAlgebra"); Pkg.add("NCDatasets"); Pkg.add("SmallCollections"); Pkg.add("Test"); Pkg.instantiate(); Pkg.status()'
```

