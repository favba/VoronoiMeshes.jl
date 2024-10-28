using VoronoiMeshes, TensorsLite, ImmutableVectors, TensorsLiteGeometry, NCDatasets, LinearAlgebra
using Test

function my_approx(a, b, atol=1e-8)
    all(x -> isapprox(x[1], x[2], atol = atol), zip(a, b))
end

fix_longitude(lon::Number) = lon > π ? lon - 2π : lon
fix_longitude(lonVec::Vector) = map(fix_longitude, lonVec)

@testset "CellInfo creation" begin
    mesh_iso = VoronoiMesh("mesh.nc")
    mesh_dist = VoronoiMesh("mesh_distorted.nc")
    mesh_spherical = VoronoiMesh("x1.4002.grid.nc")

    for mesh in (mesh_iso, mesh_dist, mesh_spherical)
        cells = mesh.cells
        if mesh === mesh_iso
            @test my_approx(cells.position, VoronoiMeshes.compute_cell_centroid(cells))
        end
        @test my_approx(cells.area, VoronoiMeshes.compute_cell_area(cells))
        @test my_approx(fix_longitude(cells.longitude), fix_longitude(VoronoiMeshes.compute_cell_longitude(cells)))
        @test my_approx(cells.latitude, VoronoiMeshes.compute_cell_latitude(cells))
    end
end

