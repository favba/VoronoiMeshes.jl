@setup_workload begin

    include("mesh_sph6_inc.jl")

    @compile_workload begin
        d = VoronoiDiagram(SphericalVoronoiDiagram(generators, vertices, verticesOnCell, cellsOnVertex, meshDensity, sphere_radius))
        m = VoronoiMesh(fix_diagram!(d))

        m.cells.area
        m.cells.centroid
        m.cells.edgesSign
        m.cells.longitude
        m.cells.latitude
        m.cells.normal
        m.cells.zonalVector
        m.cells.meridionalVector
        m.cells.areaMimetic

        m.vertices.area
        m.vertices.kiteAreas
        m.vertices.centroid
        m.vertices.edgesSign
        m.vertices.longitude
        m.vertices.latitude
        m.vertices.areaMimetic
        m.vertices.kiteAreasMimetic

        m.edges.length
        m.edges.lengthDual
        m.edges.midpoint
        m.edges.normal
        m.edges.tangent
        m.edges.angle
        m.edges.longitude
        m.edges.latitude

        graph_partition(m)

        warn_mesh_issues(m)
        
    end    
end
