abstract type VoronoiDiagram{MAX_EDGES, TI <: Integer, TF <: Real} end

struct PlanarPeriodicVoronoiDiagram{MAX_EDGES, TI, TF} <: VoronoiDiagram{MAX_EDGES, TI, TF}
    generators::Vec2DxyArray{TF, 1}
    vertices::Vec2DxyArray{TF, 1}
    verticesOnCell::ImVecArray{MAX_EDGES, TI, 1}
    cellsOnVertex::Vector{NTuple{3, TI}}
    x_period::TF
    y_period::TF
end

struct SphericalVoronoiDiagram{MAX_EDGES, TI, TF} <: VoronoiDiagram{MAX_EDGES, TI, TF}
    generators::Vec3DArray{TF, 1}
    vertices::Vec3DArray{TF, 1}
    verticesOnCell::ImVecArray{MAX_EDGES, TI, 1}
    cellsOnVertex::Vector{NTuple{3, TI}}
    sphere_radius::TF
end

