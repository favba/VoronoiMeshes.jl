module DelaunayTriangulationExt

using Zeros, VoronoiMeshes, DelaunayTriangulation, TensorsLite, ImmutableVectors, TensorsLiteGeometry

"""
    probable_dc(N::Integer, lx::Number, ly::Number, [max_densitty=1]) -> dc::Real

Given a mesh with `N` cells, (`lx`, `ly`) periods, and maximum density `max_density`, returns an educated guess at what the minimum cell distance `dc` will be.
"""
function probable_dc(N::Integer, lx::Real, ly::Real, max_den = true)
    0.95231281 * sqrt(4 * (lx * ly / N) / pi) / max_den
end

function expand_periodic_points(points::AbstractVector{<:Vec2Dxy}, xp::Number, yp::Number)
    points_indices = eachindex(points)

    N = length(points)
    #N_div = 2
    N_div = Int(floor(xp / (5 * probable_dc(N, xp, yp))))

    pxi = filter(i -> (points[i].x < xp / N_div), eachindex(points))
    px = @view(points[pxi]) .+ (xp * 𝐢)

    pxpyi = filter(i -> ((points[i].x < xp / N_div) && (points[i].y < yp / N_div)), eachindex(points))
    pxpy = @view(points[pxpyi]) .+ (xp * 𝐢 + yp * 𝐣)

    pyi = filter(i -> (points[i].y < yp / N_div), eachindex(points))
    py = @view(points[pyi]) .+ (yp * 𝐣)

    mxpyi = filter(i -> ((points[i].x > xp * (1 - 1 / N_div)) && (points[i].y < yp / N_div)), eachindex(points))
    mxpy = @view(points[mxpyi]) .+ (-xp * 𝐢 + yp * 𝐣)

    mxi = filter(i -> (points[i].x > xp * (1 - 1 / N_div)), eachindex(points))
    mx = @view(points[mxi]) .+ (-xp * 𝐢)

    mxmyi = filter(i -> ((points[i].x > xp * (1 - 1 / N_div)) && (points[i].y > yp * (1 - 1 / N_div))), eachindex(points))
    mxmy = @view(points[mxmyi]) .+ (-xp * 𝐢 - yp * 𝐣)

    myi = filter(i -> (points[i].y > yp * (1 - 1 / N_div)), eachindex(points))
    my = @view(points[myi]) .+ (-yp * 𝐣)

    pxmyi = filter(i -> ((points[i].x < xp / N_div) && (points[i].y > yp * (1 - 1 / N_div))), eachindex(points))
    pxmy = @view(points[pxmyi]) .+ (xp * 𝐢 - yp * 𝐣)

    new_points = vcat(points, px, pxpy, py, mxpy, mx, mxmy, my, pxmy)

    inds = vcat(points_indices, pxi, pxpyi, pyi, mxpyi, mxi, mxmyi, myi, pxmyi)

    return new_points, inds
end

function generate_inital_points(N::Integer, lx::Real, ly::Real)
    Nquadrant, rest = divrem(N, 9)
    if rest == 0
        Nq1 = Nq2 = Nq3 = Nq4 = Nq5 = Nq6 = Nq7 = Nq8 = Nq9 = Nquadrant
    elseif rest == 1
        Nq2 = Nq3 = Nq4 = Nq5 = Nq6 = Nq7 = Nq8 = Nq9 = Nquadrant
        Nq1 = Nq2 + 1
    elseif rest == 2
        Nq3 = Nq4 = Nq5 = Nq6 = Nq7 = Nq8 = Nq9 = Nquadrant
        Nq1 = Nq2 = Nq3 + 1
    elseif rest == 3
        Nq4 = Nq5 = Nq6 = Nq7 = Nq8 = Nq9 = Nquadrant
        Nq1 = Nq2 = Nq3 = Nq4 + 1
    elseif rest == 4
        Nq5 = Nq6 = Nq7 = Nq8 = Nq9 = Nquadrant
        Nq1 = Nq2 = Nq3 = Nq4 = Nq5 + 1
    elseif rest == 5
        Nq6 = Nq7 = Nq8 = Nq9 = Nquadrant
        Nq1 = Nq2 = Nq3 = Nq4 = Nq5 = Nq6 + 1
    elseif rest == 6
        Nq7 = Nq8 = Nq9 = Nquadrant
        Nq1 = Nq2 = Nq3 = Nq4 = Nq5 = Nq6 = Nq7 + 1
    elseif rest == 7
        Nq8 = Nq9 = Nquadrant
        Nq1 = Nq2 = Nq3 = Nq4 = Nq5 = Nq6 = Nq7 = Nq8 + 1
    elseif rest == 8
        Nq9 = Nquadrant
        Nq1 = Nq2 = Nq3 = Nq4 = Nq5 = Nq6 = Nq7 = Nq8 = Nq9 + 1
    end
    Vq1 = VecArray(x = (rand(Nq1) .* lx / 3), y = rand(Nq1) .* ly / 3)
    Vq2 = VecArray(x = (rand(Nq2) .* lx / 3 .+ lx / 3), y = rand(Nq2) .* ly / 3)
    Vq3 = VecArray(x = (rand(Nq3) .* lx / 3 .+ 2 * lx / 3), y = rand(Nq3) .* ly / 3)
    Vq4 = VecArray(x = (rand(Nq4) .* lx / 3), y = rand(Nq4) .* ly / 3 .+ ly / 3)
    Vq5 = VecArray(x = (rand(Nq5) .* lx / 3 .+ lx / 3), y = rand(Nq5) .* ly / 3 .+ ly / 3)
    Vq6 = VecArray(x = (rand(Nq6) .* lx / 3 .+ 2 * lx / 3), y = rand(Nq6) .* ly / 3 .+ ly / 3)
    Vq7 = VecArray(x = (rand(Nq7) .* lx / 3), y = rand(Nq7) .* ly / 3 .+ 2 * ly / 3)
    Vq8 = VecArray(x = (rand(Nq8) .* lx / 3 .+ lx / 3), y = rand(Nq8) .* ly / 3 .+ 2 * ly / 3)
    Vq9 = VecArray(x = (rand(Nq9) .* lx / 3 .+ 2 * lx / 3), y = rand(Nq9) .* ly / 3 .+ 2 * ly / 3)
    return vcat(Vq1, Vq2, Vq3, Vq4, Vq5, Vq6, Vq7, Vq8, Vq9)
end

@inline const_density(𝐱) = TensorsLite.Zeros.One()

function fill_with_polygon_mass_centroids!(new_points, N::Integer, lx::Number, ly::Number, voro, ρ::F = const_density) where {F <: Function}
    polygons = voro.polygons
    polygon_points = voro.polygon_points

    #Compute mass centroid for central polygons
    @parallel for i in Base.OneTo(N)
        v = polygons[i]
        vertices_position = polygon_points[@view(v[1:(end - 1)])]
        @inbounds new_points[i] = periodic_to_base_point(mass_centroid(ρ, reinterpret(Vec2Dxy{Float64}, vertices_position)), lx, ly)
    end

    return new_points
end

function centroidal_voronoi_loyd(::TT, vor::TV, initial_generator_points, N::Integer, lx::Real, ly::Real; density::F = const_density, max_iter::Integer = 20000, rtol::Real = 1.0e-4 / 100) where {TT <: Triangulation, TV <: VoronoiTessellation, F <: Function}

    initial_new = similar(initial_generator_points)

    fill_with_polygon_mass_centroids!(initial_new, N, lx, ly, vor)

    p_new, inds_new = expand_periodic_points(initial_new, lx, ly)

    new_tri::TT = triangulate(reinterpret(NTuple{2, Float64}, p_new))
    new_vor::TV = voronoi(new_tri, clip = true)

    expected_dc = probable_dc(N, lx, ly, maximum(density, initial_new))

    println("Performing Lloyd's iteration with relative tolerance = $rtol and maximum number of iterations: $max_iter")
    iter = 0
    while !(
            mapreduce(
                (x, y) -> isapprox(
                    x, y;
                    atol = rtol * expected_dc
                ),
                &, initial_generator_points, initial_new
            )
        ) &&
            iter < max_iter
        iter += 1
        mod(iter, 10) == 0 && print("Lloyd iteration number $iter \u001b[1000D")

        initial_generator_points .= initial_new
        vor::TV = new_vor

        fill_with_polygon_mass_centroids!(initial_new, N, lx, ly, vor, density)

        p_new, inds_new = expand_periodic_points(initial_new, lx, ly)

        new_tri = triangulate(reinterpret(NTuple{2, Float64}, p_new))
        new_vor = voronoi(new_tri, clip = true)
    end

    println("Stopped Lloyd's algorithm at iteration $iter")

    return initial_new, p_new, inds_new, new_vor, iter
end

is_inside_box(p::Vec2Dxy, lx::Number, ly::Number) = (0 <= p.x < lx) && (0 <= p.y < ly)

function check_and_fix_periodicity_points(points, lx::Number, ly::Number)
    x_max, _ = findmax(v -> v.x, points)
    x_min, _ = findmin(v -> v.x, points)
    y_max, _ = findmax(v -> v.y, points)
    y_min, _ = findmin(v -> v.y, points)

    if ((x_max - x_min) >= lx) || ((y_max - y_min) >= ly)
        throw(DomainError((lx, ly), "Points cannot be cast to periodic domain [0, $lx) × [0, $ly)"))
    end

    if (0 <= x_min < x_max < lx) && (0 <= y_min < y_max < ly)
        return points
    end

    displacement = (((x_max + x_min) / 2) * 𝐢 + ((y_max + y_min) / 2) * 𝐣) - ((lx / 2) * 𝐢 + (ly / 2) * 𝐣)

    @warn "Points given will be shifted by ($(displacement.x), $(displacement.y)) in order to be centered in the [0, $lx) × [0, $ly) domain"

    points_periodic = similar(points)
    points_periodic .= points .- displacement
    return points_periodic
end

function generate_periodic_centroidal_voronoi(points::AbstractVector, lx::Real, ly::Real; density::F = const_density, max_iter::Integer = 20000, rtol::Real = 1.0e-4 / 100) where {F <: Function}

    initial_generator_points = check_and_fix_periodicity_points(points, lx, ly)

    N = length(initial_generator_points)
    p, inds = expand_periodic_points(initial_generator_points, lx, ly)
    p_tuple = reinterpret(NTuple{2, Float64}, p)
    tri = triangulate(p_tuple)
    vor = voronoi(tri, clip = true)

    initial_new, p_new, inds_new, new_vor, n_iter = centroidal_voronoi_loyd(tri, vor, initial_generator_points, N, lx, ly, density = density, max_iter = max_iter, rtol = rtol)
    return new_vor
end

function find_base_point_index(point::Vec2Dxy, interior_points, x_period::Number, y_period::Number)
    findfirst(j -> isapprox_periodic(point, @inbounds(interior_points[j]), x_period, y_period), eachindex(interior_points))
end

function extract_periodic_vertices_and_cells(N::Integer, lx::Number, ly::Number, vor::VoronoiTessellation)
    #vertices position as tuples
    vertices_pos_tuple = vor.polygon_points

    #function that filters (choses) points that are inside the [0, lx) × [0, ly) domain
    fil_func = function (i)
        vp = vertices_pos_tuple[i]
        return (0 <= vp[1] < lx) && (0 <= vp[2] < ly)
    end

    #indices belonging to vertices that are in the base domain [0, lx) × [0, ly)
    interior_vertex_i = filter(fil_func, eachindex(vertices_pos_tuple))

    nVertices = length(interior_vertex_i)
    interior_vertex_pos = VecArray(
        x = map(a -> a[1], @view(vertices_pos_tuple[interior_vertex_i])),
        y = map(a -> a[2], @view(vertices_pos_tuple[interior_vertex_i]))
    )

    T = ImmutableVector{16, Int32}
    verticesOnCell_global = vor.polygons
    verticesOnCell = Vector{T}(undef, N)
    for c in eachindex(verticesOnCell)
        vocg = T(@view(verticesOnCell_global[c][1:(end - 1)]))
        vps = map(x -> reinterpret(Vec2Dxy{Float64}, getindex(vertices_pos_tuple, x)), vocg)
        verticesOnCell[c] = map(v -> find_base_point_index(v, interior_vertex_pos, lx, ly), vps)
    end

    cell_pos = parent(vor.triangulation.points)
    interior_cell_pos = cell_pos[Base.OneTo(N)]
    cellsOnVertex_global = vor.circumcenter_to_triangle
    cellsOnVertex = Vector{NTuple{3, Int32}}(undef, nVertices)
    for v in eachindex(cellsOnVertex)
        gcov = cellsOnVertex_global[interior_vertex_i[v]]
        cps = map(x -> getindex(cell_pos, x), gcov)
        cellsOnVertex[v] = map(c -> find_base_point_index(c, interior_cell_pos, lx, ly), cps)
    end

    return interior_cell_pos, interior_vertex_pos, verticesOnCell, cellsOnVertex
end

function VoronoiMeshes.PlanarVoronoiDiagram(initial_generator_points::AbstractVector{<:Vec2Dxy}, lx::Real, ly::Real; density::F = const_density, max_iter::Integer = 20000, rtol::Real = 1.0e-4 / 100) where {F <: Function}

    N = length(initial_generator_points)

    voro = generate_periodic_centroidal_voronoi(
        initial_generator_points, lx, ly;
        density = density, max_iter = max_iter, rtol = rtol
    )

    generators, vertices, verticesOnCell, cellsOnVertex = extract_periodic_vertices_and_cells(N, lx, ly, voro)

    meshDensity = Vector{nonzero_eltype(eltype(generators))}(undef, N)

    @inbounds for i in eachindex(meshDensity)
        meshDensity[i] = density(generators[i])
    end

    maxEdges = maximum(length, verticesOnCell)
    if maxEdges == 6
        verticesOnCell_6 = ImVecArray{6, Int32}(N)
        verticesOnCell_6 .= verticesOnCell
        return PlanarVoronoiDiagram(generators, vertices, verticesOnCell_6, cellsOnVertex, meshDensity, lx, ly)
    elseif maxEdges == 7
        verticesOnCell_7 = ImVecArray{7, Int32}(N)
        verticesOnCell_7 .= verticesOnCell
        return PlanarVoronoiDiagram(generators, vertices, verticesOnCell_7, cellsOnVertex, meshDensity, lx, ly)
    elseif maxEdges == 8
        verticesOnCell_8 = ImVecArray{8, Int32}(N)
        verticesOnCell_8 .= verticesOnCell
        return PlanarVoronoiDiagram(generators, vertices, verticesOnCell_8, cellsOnVertex, meshDensity, lx, ly)
    elseif maxEdges == 9
        verticesOnCell_9 = ImVecArray{9, Int32}(N)
        verticesOnCell_9 .= verticesOnCell
        return PlanarVoronoiDiagram(generators, vertices, verticesOnCell_9, cellsOnVertex, meshDensity, lx, ly)
    elseif maxEdges > 10
        throw(error("Generated Voronoi diagram has polygons of more than 10 sides"))
    else
        verticesOnCell_10 = ImVecArray{10, Int32}(N)
        verticesOnCell_10 .= verticesOnCell
        return PlanarVoronoiDiagram(generators, vertices, verticesOnCell_10, cellsOnVertex, meshDensity, lx, ly)
    end
end

function VoronoiMeshes.VoronoiDiagram(initial_generator_points::AbstractVector{<:Vec2Dxy}, lx::Real, ly::Real; density::F = const_density, max_iter::Integer = 20000, rtol::Real = 1.0e-4 / 100) where {F <: Function}
    VoronoiDiagram(PlanarVoronoiDiagram(initial_generator_points, lx, ly, density = density, max_iter = max_iter, rtol = rtol))
end

function VoronoiMeshes.PlanarVoronoiDiagram(N::Integer, lx::Real, ly::Real; density::F = const_density, max_iter::Integer = 20000, rtol::Real = 1.0e-4 / 100) where {F <: Function}
    points = generate_inital_points(N, lx, ly)
    return VoronoiMeshes.PlanarVoronoiDiagram(points, lx, ly, density = density, max_iter = max_iter, rtol = rtol)
end

function VoronoiMeshes.VoronoiDiagram(N::Integer, lx::Real, ly::Real; density::F = const_density, max_iter::Integer = 20000, rtol::Real = 1.0e-4 / 100) where {F <: Function}
    VoronoiDiagram(PlanarVoronoiDiagram(N, lx, ly, density = density, max_iter = max_iter, rtol = rtol))
end

function VoronoiMeshes.VoronoiMesh(initial_generator_points::AbstractVector{<:Vec2Dxy}, lx::Real, ly::Real; density::F = const_density, max_iter::Integer = 20000, rtol::Real = 1.0e-4 / 100) where {F <: Function}
    VoronoiMesh(VoronoiDiagram(PlanarVoronoiDiagram(initial_generator_points, lx, ly, density = density, max_iter = max_iter, rtol = rtol)))
end

function VoronoiMeshes.VoronoiMesh(N::Integer, lx::Real, ly::Real; density::F = const_density, max_iter::Integer = 20000, rtol::Real = 1.0e-4 / 100) where {F <: Function}
    VoronoiMesh(VoronoiDiagram(PlanarVoronoiDiagram(N, lx, ly, density = density, max_iter = max_iter, rtol = rtol)))
end

include("precompile_delaunay.jl")

end # Module
