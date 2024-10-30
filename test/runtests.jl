using VoronoiMeshes, TensorsLite, ImmutableVectors, TensorsLiteGeometry, NCDatasets, LinearAlgebra
using Test

function my_approx(a, b, atol = 1.0e-8)
    all(x -> isapprox(x[1], x[2], atol = atol), zip(a, b))
end

function my_approx(a::AbstractVector{<:NTuple{N}}, b::AbstractVector{<:NTuple{N}}, atol = 1.0e-8) where {N}
    all(x -> mapreduce((y, z) -> isapprox(y, z, atol = atol), &, x[1], x[2]), zip(a, b))
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
            xp = mesh.x_period
            yp = mesh.y_period
            @test my_approx(periodic_to_base_point.(cells.position, xp, yp), periodic_to_base_point.(cells.centroid, xp, yp))
        end
        if mesh === mesh_spherical
            @test my_approx(cells.position, cells.centroid, 1.0e-5)
        end
        @test my_approx(cells.area, VoronoiMeshes.compute_cell_area(cells))
        @test my_approx(fix_longitude(cells.longitude), fix_longitude(VoronoiMeshes.compute_cell_longitude(cells)))
        @test my_approx(cells.latitude, VoronoiMeshes.compute_cell_latitude(cells))

        zonal = mesh.cells.zonalVector
        meridional = mesh.cells.meridionalVector
        normal = mesh.cells.normal
        @test my_approx(normal, zonal .× meridional)
        @test my_approx(zonal, meridional .× normal)
        @test my_approx(meridional, normal .× zonal)
    end
end

@testset "VertexInfo creation" begin
    mesh_iso = VoronoiMesh("mesh.nc")
    mesh_dist = VoronoiMesh("mesh_distorted.nc")
    mesh_spherical = VoronoiMesh("x1.4002.grid.nc")

    for mesh in (mesh_iso, mesh_dist, mesh_spherical)
        vertices = mesh.vertices
        if mesh === mesh_iso
            xp = mesh.x_period
            yp = mesh.y_period
            @test my_approx(periodic_to_base_point.(vertices.position, xp, yp), periodic_to_base_point.(vertices.centroid, xp, yp))
        end
        @test my_approx(vertices.area, VoronoiMeshes.compute_vertex_area(vertices))
        @test my_approx(fix_longitude(vertices.longitude), fix_longitude(VoronoiMeshes.compute_vertex_longitude(vertices)))
        @test my_approx(vertices.latitude, VoronoiMeshes.compute_vertex_latitude(vertices))
        @test my_approx(vertices.kiteAreas, VoronoiMeshes.compute_vertex_kiteAreas(vertices))

        #zonal = mesh.vertices.zonalVector
        #meridional = mesh.vertices.meridionalVector
        #normal = mesh.vertices.normal
        #@test my_approx(normal, zonal .× meridional)
        #@test my_approx(zonal, meridional .× normal)
        #@test my_approx(meridional, normal .× zonal)
    end
end

@testset "EdgeInfo creation" begin
    mesh_iso = VoronoiMesh("mesh.nc")
    mesh_dist = VoronoiMesh("mesh_distorted.nc")
    mesh_spherical = VoronoiMesh("x1.4002.grid.nc")

    for mesh in (mesh_iso, mesh_dist, mesh_spherical)
        edges = mesh.edges
        if mesh === mesh_iso
            xp = mesh.x_period
            yp = mesh.y_period
            @test my_approx(periodic_to_base_point.(edges.position, xp, yp), periodic_to_base_point.(edges.midpoint, xp, yp))
        end
        @test my_approx(edges.length, VoronoiMeshes.compute_edge_length(edges))
        @test my_approx(edges.cellsDistance, VoronoiMeshes.compute_edge_cellsDistance(edges))
        @test my_approx(abs.(edges.angle), abs.(VoronoiMeshes.compute_edge_angle(edges)), 1.1 * π / 180)
        @test my_approx(fix_longitude(edges.longitude), fix_longitude(VoronoiMeshes.compute_edge_longitude(edges)))
        @test my_approx(edges.latitude, VoronoiMeshes.compute_edge_latitude(edges))
        @test my_approx(edges.normal, VoronoiMeshes.compute_edge_normal(edges))

        if mesh === mesh_spherical
            @test my_approx(normalize.(edges.position), edges.normal .× edges.tangent)
        else
            @test all(𝐤 .≈ edges.normal .× edges.tangent)
        end
    end
end
