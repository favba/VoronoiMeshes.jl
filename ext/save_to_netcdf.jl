using NCDatasets.CommonDataModel

function write_dim_diagram!(ds::NCDataset, diag::AbstractVoronoiDiagram{S, maxEdges}) where {S, maxEdges}
    dim = ds.dim
    dim["nCells"] = length(diag.generators)
    dim["nVertices"] = length(diag.vertices)
    dim["maxEdges"] = maxEdges
    dim["vertexDegree"] = 3
    dim["R3"] = 3
    dim["TWO"] = 2
    return ds
end

function write_diagram_fields!(ds::NCDataset, cpos, vpos, meshDensity, force3D::Bool = false)

    is3D = eltype(cpos.z) !== Zeros.Zero
    TF = eltype(cpos.x)

    defVar(ds, "xCell", cpos.x, ("nCells",); attrib = [
        "units" => "m",
        "long_name" => "Cartesian x-coordinate of cells"
    ])

    defVar(ds, "yCell", cpos.y, ("nCells",); attrib = [
        "units" => "m",
        "long_name" => "Cartesian y-coordinate of cells"
    ])

    if (is3D | force3D)
        cposz = is3D ? cpos.z : zeros(TF,length(cpos))

        defVar(ds, "zCell", cposz, ("nCells",); attrib = [
            "units" => "m",
            "long_name" => "Cartesian z-coordinate of cells"
        ])
    end

    defVar(ds, "xVertex", vpos.x, ("nVertices",); attrib = [
        "units" => "m",
        "long_name" => "Cartesian x-coordinate of vertices"
    ])
    defVar(ds, "yVertex", vpos.y, ("nVertices",); attrib = [
        "units" => "m",
        "long_name" => "Cartesian y-coordinate of vertices"
    ])

    if (is3D | force3D)
        vposz = is3D ? vpos.z : zeros(TF,length(vpos))

        defVar(ds, "zVertex", vposz, ("nVertices",); attrib = [
            "units" => "m",
            "long_name" => "Cartesian z-coordinate of vertices"
        ])
    end

    defVar(ds, "meshDensity", meshDensity, ("nCells",); attrib = [
        "units" => "-",
        "long_name" => "Mesh density function (used when generating the mesh) evaluated at a cell"
    ])

    return ds
end

function write_diagram_indexing_fields!(ds, diag::AbstractVoronoiDiagram{S, mE, TI, TF}) where {S, mE, TI, TF}
    defVar(
        ds, "cellsOnVertex", reinterpret(reshape, TI, diag.cellsOnVertex),
        ("vertexDegree", "nVertices"), attrib = 
        [
            "units" => "-",
            "long_name" => "IDs of the cells that meet at a vertex"
        ]
    )

    defVar(
        ds, "verticesOnCell", reinterpret(reshape, TI, diag.verticesOnCell.data),
        ("maxEdges", "nCells"), attrib = [
            "units" => "-",
            "long_name" => "IDs of vertices (corner points) of a cell"
        ]
    )
    defVar(
        ds, "nEdgesOnCell", TI.(diag.verticesOnCell.length),
        ("nCells",), attrib = [
            "units" => "-",
            "long_name" => "Number of edges forming the boundary of a cell"
        ]
    )
    return ds
end

function write_diagram_data!(ds::NCDataset, diag::AbstractVoronoiDiagram, force3D::Bool = false)
    write_dim_diagram!(ds, diag)
    cpos = diag.generators
    vpos = diag.vertices
    mD = diag.meshDensity
    write_diagram_fields!(ds, cpos, vpos, mD, force3D)
    write_diagram_indexing_fields!(ds, diag)
end

function save_to_netcdf!(ds::NCDataset, diag::PlanarVoronoiDiagram{maxEdges, TI}, force3D::Bool = false) where {maxEdges, TI}
    attrib = ds.attrib
    attrib["on_a_sphere"] = "NO"
    attrib["is_periodic"] = "YES"
    attrib["x_period"] = diag.x_period
    attrib["y_period"] = diag.y_period
    write_diagram_data!(ds, diag, force3D)
end

function save_to_netcdf!(ds::NCDataset, diag::SphericalVoronoiDiagram{maxEdges, TI}, force3D::Bool = false) where {maxEdges, TI}
    attrib = ds.attrib
    attrib["on_a_sphere"] = "YES"
    attrib["sphere_radius"] = diag.sphere_radius
    write_diagram_data!(ds, diag, force3D)
end

save_to_netcdf!(ds::NCDataset, diag::VoronoiDiagram, force3D::Bool = false) = save_to_netcdf!(ds, get_diagram(diag), force3D)

function save_to_netcdf(filename, diag::AbstractVoronoiDiagram; format = :netcdf4)
    NCDataset(filename, "c", format = format) do ds
        save_to_netcdf!(ds, diag)
    end
end

function write_base_cell_data!(ds::NCDataset, cells::Cells{S, mE, TI}) where {S, mE, TI}
    defVar(
        ds, "edgesOnCell", reinterpret(reshape, TI, cells.edges.data),
        ("maxEdges", "nCells"), attrib = [
            "units" => "-",
            "long_name" => "IDs of edges forming the boundary of a cell"
        ]
    )
     defVar(
        ds, "cellsOnCell", reinterpret(reshape, TI, cells.cells.data),
        ("maxEdges", "nCells"), attrib = [
            "units" => "-",
            "long_name" => "IDs of cells neighboring a cell"
        ]
    )
     return ds
end

function write_computed_cell_data!(ds::NCDataset, cells::Cells{S, mE, TI, TF}, force3D::Bool = false) where {S, mE, TI, TF}
    cinfo = cells.info
    ncells = cells.n

    if isdefined(cinfo, :centroid)
        defVar(ds, "xCellCentroid", cinfo.centroid.x, ("nCells",); attrib = [
            "units" => "m",
            "long_name" => "Cartesian x-coordinate of cells centroid"
        ])
        defVar(ds, "yCellCentroid", cinfo.centroid.y, ("nCells",); attrib = [
            "units" => "m",
            "long_name" => "Cartesian y-coordinate of cells centroid"
        ])
        if (S | force3D)
            ccz = S ? cinfo.centroid.z : zeros(TF, ncells)
            defVar(ds, "zCellCentroid", ccz, ("nCells",); attrib = [
                "units" => "m",
                "long_name" => "Cartesian z-coordinate of cells centroid"
            ])
        end
    end

    if isdefined(cinfo, :area)
        defVar(ds, "areaCell", cinfo.area, ("nCells",); attrib = [
            "units" => "m^2",
            "long_name" => ( S ? "Spherical area of a Voronoi cell" : "Area of a Voronoi cell")
        ])
    end

    if isdefined(cinfo, :longitude)
        defVar(ds, "lonCell", cinfo.longitude, ("nCells",); attrib = [
            "units" => "rad",
            "long_name" => "Longitude of cells"
        ])
    end
    if isdefined(cinfo, :latitude)
        defVar(ds, "latCell", cinfo.latitude, ("nCells",); attrib = [
            "units" => "rad",
            "long_name" => "Latitude of cells"
        ])
    end
    if isdefined(cinfo, :normal)
        cellNormal = zeros(TF, 3, ncells)
        normal = cinfo.normal
        @inbounds for c in 1:ncells
            n = normal[c]
            cellNormal[1,c] = n.x
            cellNormal[2,c] = n.y
            cellNormal[3,c] = n.z
        end
        defVar(ds, "localVerticalUnitVectors", cellNormal, ("R3", "nCells"); attrib = [
            "units" => "unitless",
            "long_name" => "Cartesian components of the vector pointing in the local vertical direction for a cell"
        ])
    end
    if (isdefined(cinfo, :zonalVector) & isdefined(cinfo, :meridionalVector))
        cellTangent = zeros(TF, 3, 2, ncells)
        zVector = cinfo.zonalVector
        mVector = cinfo.meridionalVector
        @inbounds for c in 1:ncells
            zv = zVector[c]
            cellTangent[1, 1, c] = zv.x
            cellTangent[2, 1, c] = zv.y
            cellTangent[3, 1, c] = zv.z

            mv = mVector[c]
            cellTangent[1, 2, c] = mv.x
            cellTangent[2, 2, c] = mv.y
            cellTangent[3, 2, c] = mv.z
        end
        defVar(ds, "cellTangentPlane", cellTangent, ("R3", "TWO", "nCells"); attrib = [
            "units" => "unitless",
            "long_name" => "Components of a pair of vectors defining the tangent plane at a cell"
        ])
    end
end

function write_cell_data!(ds::NCDataset, cells::Cells; force3D::Bool=false, write_computed::Bool = false)
    write_base_cell_data!(ds, cells)
    write_computed && write_computed_cell_data!(ds, cells, force3D)
    return ds
end

function write_base_vertex_data!(ds::NCDataset, vertices::Vertices{S, mE, TI}) where {S, mE, TI}
    defVar(
        ds, "edgesOnVertex", reinterpret(reshape, TI, vertices.edges),
        ("vertexDegree", "nVertices"), attrib = [
            "units" => "-",
            "long_name" => "IDs of the edges that meet at a vertex"
        ]
    )
     return ds
end

function write_computed_vertex_data!(ds::NCDataset, vertices::Vertices{S, mE, TI, TF}, force3D::Bool = false) where {S, mE, TI, TF}
    vinfo = vertices.info
    nvertices = vertices.n

    if isdefined(vinfo, :centroid)
        defVar(ds, "xVertexCentroid", vinfo.centroid.x, ("nVertices",); attrib = [
            "units" => "m",
            "long_name" => "Cartesian x-coordinate of triangles centroid"
        ])

        defVar(ds, "yVertexCentroid", vinfo.centroid.y, ("nVertices",); attrib = [
            "units" => "m",
            "long_name" => "Cartesian y-coordinate of triangles centroid"
        ])
        if (S | force3D)
            cz = S ? vinfo.centroid.z : zeros(TF, nvertices)
            defVar(ds, "zVertexCentroid", cz, ("nVertices",); attrib = [
                "units" => "m",
                "long_name" => "Cartesian z-coordinate of triangles centroid"
            ])
        end
    end

    if isdefined(vinfo, :area)
        defVar(ds, "areaTriangle", vinfo.area, ("nVertices",); attrib = [
            "units" => "m^2",
            "long_name" => ( S ? "Spherical area of Delaunay Triangle" : "Area of a Delaunay Triangle")
        ])
    end

    if isdefined(vinfo, :kiteAreas)
        defVar(ds, "kiteAreasOnVertex", reinterpret(reshape, TF, vinfo.kiteAreas), ("vertexDegree", "nVertices"); attrib = [
            "units" => "m^2",
            "long_name" => "Intersection areas between primal (Voronoi) and dual (triangular) mesh cells"
        ])
    end

    if isdefined(vinfo, :longitude)
        defVar(ds, "lonVertex", vinfo.longitude, ("nVertices",); attrib = [
            "units" => "rad",
            "long_name" => "Longitude of vertices"
        ])
    end
    if isdefined(vinfo, :latitude)
        defVar(ds, "latVertex", vinfo.latitude, ("nVertices",); attrib = [
            "units" => "rad",
            "long_name" => "Latitude of vertices"
        ])
    end

end

function write_vertex_data!(ds::NCDataset, vertices::Vertices; force3D::Bool=false, write_computed::Bool = false)
    write_base_vertex_data!(ds, vertices)
    write_computed && write_computed_vertex_data!(ds, vertices, force3D)
    return ds
end

function write_base_edge_data!(ds::NCDataset, edges::Edges{S, mE, TI, TF}, force3D::Bool = false) where {S, mE, TI, TF}
    defVar(
        ds, "verticesOnEdge", reinterpret(reshape, TI, edges.vertices),
        ("TWO", "nEdges"), attrib = [
            "units" => "-",
            "long_name" => "IDs of the two vertex endpoints of an edge"
        ]
    )

    defVar(
        ds, "cellsOnEdge", reinterpret(reshape, TI, edges.cells),
        ("TWO", "nEdges"), attrib = [
            "units" => "-",
            "long_name" => "IDs of cells divided by an edge"
        ]
    )

    defVar(ds, "xEdge", edges.position.x, ("nEdges",); attrib = [
        "units" => "m",
        "long_name" => "Cartesian x-coordinate of edges"
    ])

    defVar(ds, "yEdge", edges.position.y, ("nEdges",); attrib = [
        "units" => "m",
        "long_name" => "Cartesian y-coordinate of edges"
    ])

    if (S | force3D)
        eposz = S ? edges.position.z : zeros(TF,length(edges.position))

        defVar(ds, "zEdge", eposz, ("nEdges",); attrib = [
            "units" => "m",
            "long_name" => "Cartesian z-coordinate of edges"
        ])
    end

     return ds
end

function write_computed_edge_data!(ds::NCDataset, edges::Edges{S, mE, TI, TF}, force3D::Bool = false) where {S, mE, TI, TF}

    einfo = edges.info
    nEdges = edges.n

    if isdefined(einfo, :midpoint)
        defVar(ds, "xEdgeMidpoint", einfo.midpoint.x, ("nEdges",); attrib = [
            "units" => "m",
            "long_name" => "Cartesian x-coordinate of edges midpoint"
        ])

        defVar(ds, "yEdgeMidpoint", einfo.midpoint.y, ("nEdges",); attrib = [
            "units" => "m",
            "long_name" => "Cartesian y-coordinate of edges midpoint"
        ])

        if (S | force3D)
            eposz = S ? einfo.midpoint.z : zeros(TF, nEdges)

            defVar(ds, "zEdgeMidpoint", eposz, ("nEdges",); attrib = [
                "units" => "m",
                "long_name" => "Cartesian z-coordinate of edges midpoint"
            ])
        end
    end

    if isdefined(einfo, :length)
        defVar(ds, "dvEdge", einfo.length, ("nEdges",); attrib = [
            "units" => "m",
            "long_name" => (S ? "Spherical distance between vertex endpoints of an edge" : "Distance between vertex endpoints of an edge")
        ])
    end

    if isdefined(einfo, :cellsDistance)
        defVar(ds, "dcEdge", einfo.cellsDistance, ("nEdges",); attrib = [
            "units" => "m",
            "long_name" => (S ? "Spherical distance between cells separated by an edge" : "Distance between cells separated by an edge")
        ])
    end

    if isdefined(einfo, :angle)
        defVar(ds, "angleEdge", einfo.angle, ("nEdges",); attrib = [
            "units" => "rad",
            "long_name" => "Angle between local north and the positive tangential direction of an edge"
        ])
    end

    if isdefined(einfo, :longitude)
        defVar(ds, "lonEdge", einfo.longitude, ("nEdges",); attrib = [
            "units" => "rad",
            "long_name" => "Longitude of edges"
        ])
    end
    if isdefined(einfo, :latitude)
        defVar(ds, "latEdge", einfo.latitude, ("nEdges",); attrib = [
            "units" => "rad",
            "long_name" => "Latitude of edges"
        ])
    end

    if isdefined(einfo, :normal)
        edgeNormal = zeros(TF, 3, nEdges)
        normal = einfo.normal
        @inbounds for e in 1:nEdges
            n = normal[e]
            edgeNormal[1,e] = n.x
            edgeNormal[2,e] = n.y
            edgeNormal[3,e] = n.z
        end
        defVar(ds, "edgeNormalVectors", edgeNormal, ("R3", "nEdges"); attrib = [
            "units" => "unitless",
            "long_name" => "Cartesian components of the vector normal to an edge and tangential to the surface of the sphere"
        ])
    end

    if isdefined(einfo, :tangent)
        edgeTangent = zeros(TF, 3, nEdges)
        tangent = einfo.tangent
        @inbounds for e in 1:nEdges
            t = tangent[e]
            edgeTangent[1,e] = t.x
            edgeTangent[2,e] = t.y
            edgeTangent[3,e] = t.z
        end
        defVar(ds, "edgeTangentialVectors", edgeTangent, ("R3", "nEdges"); attrib = [
            "units" => "unitless",
            "long_name" => "Cartesian components of the vector tangential to an edge and tangential to the surface of the sphere"
        ])
    end

    return ds
end

function write_edge_data!(ds::NCDataset, edges::Edges; force3D::Bool = false, write_computed::Bool = false)
    write_base_edge_data!(ds, edges, force3D)
    write_computed && write_computed_edge_data!(ds, edges, force3D)
    return ds
end

function save_to_netcdf!(ds::NCDataset, mesh::VoronoiMesh; force3D::Bool = false, write_computed::Bool= false)
    save_to_netcdf!(ds, mesh.diagram, force3D)
    write_cell_data!(ds, mesh.cells; force3D = force3D, write_computed = write_computed)
    write_vertex_data!(ds, mesh.vertices; force3D = force3D, write_computed = write_computed)
    write_edge_data!(ds, mesh.edges; force3D = force3D, write_computed = write_computed)
end

function save_to_netcdf(filename, mesh::VoronoiMesh; force3D::Bool=false, write_computed::Bool=false, format = :netcdf4)
    NCDataset(filename, "c", format = format) do ds
        save_to_netcdf!(ds, mesh; force3D = force3D, write_computed = write_computed)
    end
end

