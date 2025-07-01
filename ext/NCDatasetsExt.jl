module NCDatasetsExt

using VoronoiMeshes, NCDatasets, TensorsLite, TensorsLite.Zeros, SmallCollections
import VoronoiMeshes: save_to_netcdf, save_to_netcdf!, copy_matrix_to_fixedvector_vector!
using PrecompileTools

function on_a_sphere(ncfile::NCDatasets.NCDataset)
    oas = lowercase(strip(ncfile.attrib["on_a_sphere"]::String))
    return oas in ("yes", "y")
end

precompile(on_a_sphere, (NCDatasets.NCDataset{Nothing, Missing},))

function VoronoiMeshes.PlanarVoronoiDiagram(::Val{N}, nEdges::Vector{Int16}, ncfile::NCDatasets.NCDataset) where {N}
    verticesOnCellArray = (ncfile["verticesOnCell"][:, :])::Matrix{Int32}
    nCells = length(nEdges)

    verticesOnCell = SmallVectorArray(Vector{FixedVector{N, Int32}}(undef, nCells), nEdges)
    t1 = Threads.@spawn copy_matrix_to_fixedvector_vector!($(verticesOnCell.data), $verticesOnCellArray)

    cellsOnVertexArray = ncfile["cellsOnVertex"][:, :]::Matrix{Int32}
    nVertices = size(cellsOnVertexArray, 2)
    cellsOnVertex = Vector{FixedVector{3, Int32}}(undef, nVertices)
    t2 = Threads.@spawn copy_matrix_to_fixedvector_vector!($cellsOnVertex, $cellsOnVertexArray)

    xCell = ncfile["xCell"][:]::Vector{Float64}
    yCell = ncfile["yCell"][:]::Vector{Float64}
    generators = VecArray(x = xCell, y = yCell)

    xVertex = ncfile["xVertex"][:]::Vector{Float64}
    yVertex = ncfile["yVertex"][:]::Vector{Float64}
    vertices = VecArray(x = xVertex, y = yVertex)

    x_period = ncfile.attrib["x_period"]::Float64
    y_period = ncfile.attrib["y_period"]::Float64

    meshDensity = ncfile["meshDensity"][:]::Vector{Float64}

    wait(t1)
    wait(t2)

    return PlanarVoronoiDiagram(generators, vertices, verticesOnCell, cellsOnVertex, meshDensity, x_period, y_period)
end

for N in 6:10
    precompile(VoronoiMeshes.PlanarVoronoiDiagram, (Val{N}, Vector{Int16}, NCDatasets.NCDataset{Nothing, Missing}))
end

function VoronoiMeshes.PlanarVoronoiDiagram(ncfile::NCDatasets.NCDataset)
    on_a_sphere(ncfile) && throw(error("Mesh is not planar"))

    nEdgesOnCell = Int16.(ncfile["nEdgesOnCell"][:]::Vector{Int32})
    maxEdges = Int(maximum(nEdgesOnCell))

    #Avoid dynamic dispatch for most common cases
    if maxEdges == 6
        return PlanarVoronoiDiagram(Val{6}(), nEdgesOnCell, ncfile)
    elseif maxEdges == 7
        return PlanarVoronoiDiagram(Val{7}(), nEdgesOnCell, ncfile)
    elseif maxEdges == 8
        return PlanarVoronoiDiagram(Val{8}(), nEdgesOnCell, ncfile)
    elseif maxEdges == 9
        return PlanarVoronoiDiagram(Val{9}(), nEdgesOnCell, ncfile)
    elseif maxEdges == 10
        return PlanarVoronoiDiagram(Val{10}(), nEdgesOnCell, ncfile)
    else
        return PlanarVoronoiDiagram(Val(maxEdges), nEdgesOnCell, ncfile)
    end

end

function VoronoiMeshes.PlanarVoronoiDiagram(n::Val{maxEdges}, ncfile::NCDatasets.NCDataset) where {maxEdges}
    on_a_sphere(ncfile) && throw(error("Mesh is not planar"))

    nEdgesOnCell = Int16.(ncfile["nEdgesOnCell"][:]::Vector{Int32})
    maxEdges == Int(maximum(nEdgesOnCell)) || throw(error("nEdges not consistent with data in NetCDF file"))

    return PlanarVoronoiDiagram(n, nEdgesOnCell, ncfile)
end

for N in 6:10
    precompile(VoronoiMeshes.PlanarVoronoiDiagram, (Val{N}, NCDatasets.NCDataset{Nothing, Missing}))
end

function VoronoiMeshes.SphericalVoronoiDiagram(::Val{N}, nEdges::Vector{Int16}, ncfile::NCDatasets.NCDataset) where {N}
    verticesOnCellArray = (ncfile["verticesOnCell"][:, :])::Matrix{Int32}
    nCells = length(nEdges)

    verticesOnCell = SmallVectorArray(Vector{FixedVector{N, Int32}}(undef, nCells), nEdges)
    t1 = Threads.@spawn copy_matrix_to_fixedvector_vector!($(verticesOnCell.data), $verticesOnCellArray)

    cellsOnVertexArray = ncfile["cellsOnVertex"][:, :]::Matrix{Int32}
    nVertices = size(cellsOnVertexArray, 2)
    cellsOnVertex = Vector{FixedVector{3, Int32}}(undef, nVertices)
    t2 = Threads.@spawn copy_matrix_to_fixedvector_vector!($cellsOnVertex, $cellsOnVertexArray)

    xCell = ncfile["xCell"][:]::Vector{Float64}
    yCell = ncfile["yCell"][:]::Vector{Float64}
    zCell = ncfile["zCell"][:]::Vector{Float64}
    generators = VecArray(xCell, yCell, zCell)

    xVertex = ncfile["xVertex"][:]::Vector{Float64}
    yVertex = ncfile["yVertex"][:]::Vector{Float64}
    zVertex = ncfile["zVertex"][:]::Vector{Float64}
    vertices = VecArray(xVertex, yVertex, zVertex)

    sphere_radius = ncfile.attrib["sphere_radius"]::Float64

    meshDensity = ncfile["meshDensity"][:]::Vector{Float64}

    wait(t1)
    wait(t2)

    return SphericalVoronoiDiagram(generators, vertices, verticesOnCell, cellsOnVertex, meshDensity, sphere_radius)
end

for N in 6:10
    precompile(VoronoiMeshes.SphericalVoronoiDiagram, (Val{N}, Vector{Int16}, NCDatasets.NCDataset{Nothing, Missing}))
end

function VoronoiMeshes.SphericalVoronoiDiagram(ncfile::NCDatasets.NCDataset)
    on_a_sphere(ncfile) || throw(error("Mesh is not spherical"))

    nEdgesOnCell = Int16.(ncfile["nEdgesOnCell"][:]::Vector{Int32})
    maxEdges = Int(maximum(nEdgesOnCell))

    #Avoid dynamic dispatch for most common cases
    if maxEdges == 6
        return SphericalVoronoiDiagram(Val{6}(), nEdgesOnCell, ncfile)
    elseif maxEdges == 7
        return SphericalVoronoiDiagram(Val{7}(), nEdgesOnCell, ncfile)
    elseif maxEdges == 8
        return SphericalVoronoiDiagram(Val{8}(), nEdgesOnCell, ncfile)
    elseif maxEdges == 9
        return SphericalVoronoiDiagram(Val{9}(), nEdgesOnCell, ncfile)
    elseif maxEdges == 10
        return SphericalVoronoiDiagram(Val{10}(), nEdgesOnCell, ncfile)
    else
        return SphericalVoronoiDiagram(Val(maxEdges), nEdgesOnCell, ncfile)
    end

end

function VoronoiMeshes.SphericalVoronoiDiagram(n::Val{maxEdges}, ncfile::NCDatasets.NCDataset) where {maxEdges}
    on_a_sphere(ncfile) || throw(error("Mesh is not spherical"))

    nEdgesOnCell = Int16.(ncfile["nEdgesOnCell"][:]::Vector{Int32})
    maxEdges == Int(maximum(nEdgesOnCell)) || throw(error("nEdges not consistent with data in NetCDF file"))

    return SphericalVoronoiDiagram(n, nEdgesOnCell, ncfile)
end

for N in 6:10
    precompile(VoronoiMeshes.SphericalVoronoiDiagram, (Val{N}, NCDatasets.NCDataset{Nothing, Missing}))
end

function VoronoiMeshes.VoronoiDiagram(ncfile::NCDatasets.NCDataset)
    if on_a_sphere(ncfile)
        return VoronoiDiagram(SphericalVoronoiDiagram(ncfile))
    else
        return VoronoiDiagram(PlanarVoronoiDiagram(ncfile))
    end
end

function VoronoiMeshes.VoronoiDiagram(n::Val{maxEdges}, ncfile::NCDatasets.NCDataset) where {maxEdges}
    nEdges = Int16.(ncfile["nEdgesOnCell"][:]::Vector{Int32})
    maxEdges == Int(maximum(nEdges)) || throw(error("nEdges not consistent with data in NetCDF file"))
    if on_a_sphere(ncfile)
        return VoronoiDiagram(SphericalVoronoiDiagram(n, nEdges, ncfile))
    else
        return VoronoiDiagram(PlanarVoronoiDiagram(n, nEdges, ncfile))
    end
end

function VoronoiMeshes.CellInfo(voro::VoronoiDiagram{S, NE, TI, TF, Tz}, ncfile::NCDatasets.NCDataset) where {S, NE, TI, TF, Tz}

    cell_info = CellInfo(voro)

    haskey(ncfile, "centroidCell") && (cell_info.area = ncfile["centroidCell"][:]::Vector{Float64})
    haskey(ncfile, "areaCell") && (cell_info.area = ncfile["areaCell"][:]::Vector{Float64})
    haskey(ncfile, "lonCell") && (cell_info.longitude = ncfile["lonCell"][:]::Vector{Float64})
    haskey(ncfile, "latCell") && (cell_info.latitude = ncfile["latCell"][:]::Vector{Float64})

    if haskey(ncfile, "localVerticalUnitVectors")
        lvuva = ncfile["localVerticalUnitVectors"]
        n_z = lvuva[3, :]
        if S
            n_x = lvuva[1, :]
            n_y = lvuva[2, :]
            cell_info.normal = VecArray(n_x, n_y, n_z)
        else
            cell_info.normal = VecArray(z = n_z)
        end
    end

    if haskey(ncfile, "cellTangentPlane")
        ctp = ncfile["cellTangentPlane"]
        tux = ctp[1, 1, :]
        tvy = ctp[2, 2, :]
        if S
            tuy = ctp[2, 1, :]
            tvx = ctp[1, 2, :]
            tvz = ctp[3, 2, :]
            cell_info.zonalVector = VecArray(x = tux, y = tuy)
            cell_info.meridionalVector = VecArray(tvx, tvy, tvz)
        else
            cell_info.zonalVector = VecArray(x = tux)
            cell_info.meridionalVector = VecArray(y = tvy)
        end
    end

    return cell_info
end

VoronoiMeshes.CellInfo(ncfile::NCDatasets.NCDataset) = CellInfo(VoronoiDiagram(ncfile), ncfile)
VoronoiMeshes.CellInfo(v::Val{N}, ncfile::NCDatasets.NCDataset) where {N} = CellInfo(VoronoiDiagram(v, ncfile), ncfile)


function VoronoiMeshes.Cells(voro::VoronoiDiagram{S, NE, TI, TF, Tz}, ncfile::NCDatasets.NCDataset) where {S, NE, TI, TF, Tz}
    position = voro.generators
    n = length(position)
    vertices = voro.verticesOnCell
    nEdges = vertices.length
    edges = SmallVectorArray(similar(vertices.data), nEdges)
    copy_matrix_to_fixedvector_vector!(edges.data, ncfile["edgesOnCell"][:, :]::Matrix{Int32})
    cells = SmallVectorArray(similar(vertices.data), nEdges)
    copy_matrix_to_fixedvector_vector!(cells.data, ncfile["cellsOnCell"][:, :]::Matrix{Int32})

    return Cells(n, position, nEdges, vertices, edges, cells, CellInfo(voro, ncfile))
end

VoronoiMeshes.Cells(ncfile::NCDatasets.NCDataset) = Cells(VoronoiDiagram(ncfile), ncfile)
VoronoiMeshes.Cells(v::Val{N}, ncfile::NCDatasets.NCDataset) where {N} = Cells(VoronoiDiagram(v, ncfile), ncfile)

function VoronoiMeshes.VertexInfo(voro::VoronoiDiagram{S, NE, TI, TF, Tz}, ncfile::NCDatasets.NCDataset) where {S, NE, TI, TF, Tz}
    vertex_info = VertexInfo(voro)

    haskey(ncfile, "centroidVertex") && (vertex_info.centroid = ncfile["centroidVertex"][:]::Vector{Float64})
    haskey(ncfile, "centroidTriangle") && (vertex_info.centroid = ncfile["centroidTriangle"][:]::Vector{Float64})
    haskey(ncfile, "areaVertex") && (vertex_info.area = ncfile["areaVertex"][:]::Vector{Float64})
    haskey(ncfile, "areaTriangle") && (vertex_info.area = ncfile["areaTriangle"][:]::Vector{Float64})
    haskey(ncfile, "lonVertex") && (vertex_info.longitude = ncfile["lonVertex"][:]::Vector{Float64})
    haskey(ncfile, "latVertex") && (vertex_info.latitude = ncfile["latVertex"][:]::Vector{Float64})

    if haskey(ncfile, "kiteAreasOnVertex")
        kiteAreasArray = ncfile["kiteAreasOnVertex"][:, :]
        kiteAreas = [(kiteAreasArray[1, k], kiteAreasArray[2, k], kiteAreasArray[3, k]) for k in axes(kiteAreasArray, 2)]
        vertex_info.kiteAreas = kiteAreas
    end

    return vertex_info
end

VoronoiMeshes.VertexInfo(ncfile::NCDatasets.NCDataset) = VertexInfo(VoronoiDiagram(ncfile), ncfile)
VoronoiMeshes.VertexInfo(v::Val{N}, ncfile::NCDatasets.NCDataset) where {N} = VertexInfo(VoronoiDiagram(v, ncfile), ncfile)

function VoronoiMeshes.Vertices(voro::VoronoiDiagram{S, NE, TI, TF, Tz}, ncfile::NCDatasets.NCDataset) where {S, NE, TI, TF, Tz}
    position = voro.vertices
    n = length(position)
    cells = voro.cellsOnVertex
    edges = similar(cells)
    copy_matrix_to_fixedvector_vector!(edges, ncfile["edgesOnVertex"][:, :]::Matrix{Int32})
    return Vertices(n, position, edges, cells, VertexInfo(voro, ncfile))
end

VoronoiMeshes.Vertices(ncfile::NCDatasets.NCDataset) = Vertices(VoronoiDiagram(ncfile), ncfile)
VoronoiMeshes.Vertices(v::Val{N}, ncfile::NCDatasets.NCDataset) where {N} = Vertices(VoronoiDiagram(v, ncfile), ncfile)

function VoronoiMeshes.EdgeInfo(voro::VoronoiDiagram{S, NE, TI, TF, Tz}, ncfile::NCDatasets.NCDataset) where {S, NE, TI, TF, Tz}
    edge_info = EdgeInfo(voro)

    haskey(ncfile, "midpointEdge") && (edge_info.midpoint = ncfile["midpointEdge"][:]::Vector{Float64})
    haskey(ncfile, "lengthEdge") && (edge_info.length = ncfile["lengthEdge"][:]::Vector{Float64})
    haskey(ncfile, "dvEdge") && (edge_info.length = ncfile["dvEdge"][:]::Vector{Float64})
    haskey(ncfile, "lengthDualEdge") && (edge_info.lengthDual = ncfile["lengthDualEdge"][:]::Vector{Float64})
    haskey(ncfile, "dcEdge") && (edge_info.lengthDual = ncfile["dcEdge"][:]::Vector{Float64})
    haskey(ncfile, "angleEdge") && (edge_info.angle = ncfile["angleEdge"][:]::Vector{Float64})
    haskey(ncfile, "lonEdge") && (edge_info.longitude = ncfile["lonEdge"][:]::Vector{Float64})
    haskey(ncfile, "latEdge") && (edge_info.latitude = ncfile["latEdge"][:]::Vector{Float64})

    if haskey(ncfile, "edgeNormalVectors")
        env = ncfile["edgeNormalVectors"]
        if S
            edge_info.normal = VecArray(x = env[1, :], y = env[2, :], z = env[3, :])
        else
            edge_info.normal = VecArray(x = env[1, :], y = env[2, :])
        end
    end

    if haskey(ncfile, "edgeTangentialVectors")
        etv = ncfile["edgeTangentialVectors"]
        if S
            edge_info.tangent = VecArray(x = etv[1, :], y = etv[2, :], z = etv[3, :])
        else
            edge_info.tangent = VecArray(x = etv[1, :], y = etv[2, :])
        end
    end

    return edge_info
end

VoronoiMeshes.EdgeInfo(ncfile::NCDatasets.NCDataset) = EdgeInfo(VoronoiDiagram(ncfile), ncfile)
VoronoiMeshes.EdgeInfo(v::Val{N}, ncfile::NCDatasets.NCDataset) where {N} = EdgeInfo(VoronoiDiagram(v, ncfile), ncfile)

function VoronoiMeshes.Edges(voro::VoronoiDiagram{S, NE, TI, TF, Tz}, ncfile::NCDatasets.NCDataset) where {S, NE, TI, TF, Tz}

    verticesOnEdgeArray = ncfile["verticesOnEdge"][:, :]::Matrix{Int32}
    n = size(verticesOnEdgeArray, 2)
    vertices = Vector{FixedVector{2, Int32}}(undef, n)
    copy_matrix_to_fixedvector_vector!(vertices, verticesOnEdgeArray)

    cellsOnEdgeArray = ncfile["cellsOnEdge"][:, :]::Matrix{Int32}
    cells = Vector{FixedVector{2, Int32}}(undef, n)
    copy_matrix_to_fixedvector_vector!(cells, cellsOnEdgeArray)

    xEdge = ncfile["xEdge"][:]::Vector{Float64}
    yEdge = ncfile["yEdge"][:]::Vector{Float64}

    if S
        zEdge = ncfile["zEdge"][:]::Vector{Float64}
        position = VecArray(xEdge, yEdge, zEdge)
    else
        position = VecArray(x = xEdge, y = yEdge)
    end

    return Edges(n, position, vertices, cells, EdgeInfo(voro, ncfile))
end

VoronoiMeshes.Edges(ncfile::NCDatasets.NCDataset) = Edges(VoronoiDiagram(ncfile), ncfile)
VoronoiMeshes.Edges(v::Val{N}, ncfile::NCDatasets.NCDataset) where {N} = Edges(VoronoiDiagram(v, ncfile), ncfile)

for S in (true, false)
    for NE in 6:10
        for TI in (Int32, Int64)
            for TF in (Float32, Float64)
                for Tz in (TF, Zero)
                    precompile(VoronoiMeshes.CellInfo, (VoronoiDiagram{S, NE, TI, TF, Tz}, NCDatasets.NCDataset{Nothing, Missing}))
                    precompile(VoronoiMeshes.Cells, (VoronoiDiagram{S, NE, TI, TF, Tz}, NCDatasets.NCDataset{Nothing, Missing}))
                    precompile(VoronoiMeshes.VertexInfo, (VoronoiDiagram{S, NE, TI, TF, Tz}, NCDatasets.NCDataset{Nothing, Missing}))
                    precompile(VoronoiMeshes.Vertices, (VoronoiDiagram{S, NE, TI, TF, Tz}, NCDatasets.NCDataset{Nothing, Missing}))
                    precompile(VoronoiMeshes.EdgeInfo, (VoronoiDiagram{S, NE, TI, TF, Tz}, NCDatasets.NCDataset{Nothing, Missing}))
                    precompile(VoronoiMeshes.Edges, (VoronoiDiagram{S, NE, TI, TF, Tz}, NCDatasets.NCDataset{Nothing, Missing}))
                end
            end
        end
    end
end

function _VoronoiMesh(ncfile::NCDatasets.NCDataset)
    maxEdges = Int(maximum(ncfile["nEdgesOnCell"][:]::Vector{Int32}))
    S = on_a_sphere(ncfile)
    if S
        if maxEdges == 6
            diag6 = VoronoiDiagram(SphericalVoronoiDiagram(Val{6}(), ncfile))
            return VoronoiMesh(Cells(diag6, ncfile), Vertices(diag6, ncfile), Edges(diag6, ncfile))
        elseif maxEdges == 7
            diag7 = VoronoiDiagram(SphericalVoronoiDiagram(Val{7}(), ncfile))
            return VoronoiMesh(Cells(diag7, ncfile), Vertices(diag7, ncfile), Edges(diag7, ncfile))
        elseif maxEdges == 8
            diag8 = VoronoiDiagram(SphericalVoronoiDiagram(Val{8}(), ncfile))
            return VoronoiMesh(Cells(diag8, ncfile), Vertices(diag8, ncfile), Edges(diag8, ncfile))
        elseif maxEdges == 9
            diag9 = VoronoiDiagram(SphericalVoronoiDiagram(Val{9}(), ncfile))
            return VoronoiMesh(Cells(diag9, ncfile), Vertices(diag9, ncfile), Edges(diag9, ncfile))

        elseif maxEdges == 10
            diag10 = VoronoiDiagram(SphericalVoronoiDiagram(Val{10}(), ncfile))
            return VoronoiMesh(Cells(diag10, ncfile), Vertices(diag10, ncfile), Edges(diag10, ncfile))
        else
            diag = VoronoiDiagram(SphericalVoronoiDiagram(ncfile))
            return VoronoiMesh(Cells(diag, ncfile), Vertices(diag, ncfile), Edges(diag, ncfile))::(VoronoiMesh{true, N, Int32, Float64, Float64} where {N})
        end
    else
        if maxEdges == 6
            pdiag6 = VoronoiDiagram(PlanarVoronoiDiagram(Val{6}(), ncfile))
            return VoronoiMesh(Cells(pdiag6, ncfile), Vertices(pdiag6, ncfile), Edges(pdiag6, ncfile))
        elseif maxEdges == 7
            pdiag7 = VoronoiDiagram(PlanarVoronoiDiagram(Val{7}(), ncfile))
            return VoronoiMesh(Cells(pdiag7, ncfile), Vertices(pdiag7, ncfile), Edges(pdiag7, ncfile))
        elseif maxEdges == 8
            pdiag8 = VoronoiDiagram(PlanarVoronoiDiagram(Val{8}(), ncfile))
            return VoronoiMesh(Cells(pdiag8, ncfile), Vertices(pdiag8, ncfile), Edges(pdiag8, ncfile))
        elseif maxEdges == 9
            pdiag9 = VoronoiDiagram(PlanarVoronoiDiagram(Val{9}(), ncfile))
            return VoronoiMesh(Cells(pdiag9, ncfile), Vertices(pdiag9, ncfile), Edges(pdiag9, ncfile))

        elseif maxEdges == 10
            pdiag10 = VoronoiDiagram(PlanarVoronoiDiagram(Val{10}(), ncfile))
            return VoronoiMesh(Cells(pdiag10, ncfile), Vertices(pdiag10, ncfile), Edges(pdiag10, ncfile))
        else
            pdiag = VoronoiDiagram(PlanarVoronoiDiagram(ncfile))
            return VoronoiMesh(Cells(pdiag, ncfile), Vertices(pdiag, ncfile), Edges(pdiag, ncfile))::(VoronoiMesh{false, N, Int32, Float64, Zeros.Zero} where {N})
        end
    end
end

function VoronoiMeshes.VoronoiMesh(ncfile::NCDatasets.NCDataset, warn_issues::Bool=true)
    mesh = _VoronoiMesh(ncfile)
    if warn_issues
        Threads.@spawn VoronoiMeshes.warn_mesh_issues($mesh)
    end
    return mesh
end

function VoronoiMeshes.VoronoiMesh(v::Val{NE}, ncfile::NCDatasets.NCDataset, warn_issues::Bool = true) where {NE}
    diag = VoronoiDiagram(v, ncfile)
    mesh = VoronoiMesh(Cells(diag, ncfile), Vertices(diag, ncfile), Edges(diag, ncfile))
    if warn_issues
        Threads.@spawn VoronoiMeshes.warn_mesh_issues($mesh)
    end
    return mesh
end

for func in (
        :PlanarVoronoiDiagram, :SphericalVoronoiDiagram, :VoronoiDiagram,
        :CellInfo, :Cells, :VertexInfo, :Vertices, :EdgeInfo, :Edges,
       # :VoronoiMesh,
    )
    @eval begin
        function VoronoiMeshes.$func(file_name::String)
            f = NCDataset(file_name)
            try
                VoronoiMeshes.$func(f)
            finally
                close(f)
            end
        end
        function VoronoiMeshes.$func(v::Val, file_name::String)
            f = NCDataset(file_name)
            try
                VoronoiMeshes.$func(v, f)
            finally
                close(f)
            end
        end
        precompile(VoronoiMeshes.$func, (NCDatasets.NCDataset{Nothing, Missing},))
        precompile(VoronoiMeshes.$func, (String,))
    end
    for N in 6:10
        @eval precompile(VoronoiMeshes.$func, (Val{$N}, NCDatasets.NCDataset{Nothing, Missing}))
        @eval precompile(VoronoiMeshes.$func, (Val{$N}, String))
    end
end

function VoronoiMeshes.VoronoiMesh(file_name::String, warn_issues::Bool = true)
    f = NCDataset(file_name)
    try
        VoronoiMeshes.VoronoiMesh(f, warn_issues)
    finally
        close(f)
    end
end
function VoronoiMeshes.VoronoiMesh(v::Val, file_name::String, warn_issues::Bool = true)
    f = NCDataset(file_name)
    try
        VoronoiMeshes.VoronoiMesh(v, f, warn_issues)
    finally
        close(f)
    end
end

precompile(VoronoiMeshes.VoronoiMesh, (NCDatasets.NCDataset{Nothing, Missing}, Bool))
precompile(VoronoiMeshes.VoronoiMesh, (String, Bool))
precompile(VoronoiMeshes.VoronoiMesh, (NCDatasets.NCDataset{Nothing, Missing},))
precompile(VoronoiMeshes.VoronoiMesh, (String,))

for N in 6:10
    @eval precompile(VoronoiMeshes.VoronoiMesh, (Val{$N}, NCDatasets.NCDataset{Nothing, Missing}, Bool))
    @eval precompile(VoronoiMeshes.VoronoiMesh, (Val{$N}, String, Bool))
    @eval precompile(VoronoiMeshes.VoronoiMesh, (Val{$N}, NCDatasets.NCDataset{Nothing, Missing}))
    @eval precompile(VoronoiMeshes.VoronoiMesh, (Val{$N}, String))
end

include("save_to_netcdf.jl")

@compile_workload begin
    m_iso = VoronoiMesh("../test/mesh.nc")
    m_sphe = VoronoiMesh("../test/spherical_grid_500km.nc")
    m_dist = VoronoiMesh(fix_diagram!(VoronoiDiagram("../test/mesh_distorted_issues.nc")))
    m_sphe_2 = VoronoiMesh(fix_diagram!(VoronoiDiagram("../test/x1.4002.grid.nc")))

    save("test_save.nc", m_iso; force3D=true, write_computed=true)
    Base.Filesystem.rm("test_save.nc")
    Base.Filesystem.rm("test_save.graph.info")
    save("test_save.nc", m_sphe; write_computed=true)
    Base.Filesystem.rm("test_save.nc")
    Base.Filesystem.rm("test_save.graph.info")
end

end # module
