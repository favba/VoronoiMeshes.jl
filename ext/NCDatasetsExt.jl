module NCDatasetsExt

using VoronoiMeshes, NCDatasets, TensorsLite, TensorsLite.Zeros, ImmutableVectors

function on_a_sphere(ncfile::NCDatasets.NCDataset)
    oas = lowercase(strip(ncfile.attrib["on_a_sphere"]::String))
    return oas in ("yes","y")
end

precompile(on_a_sphere, (NCDatasets.NCDataset{Nothing, Missing},))

function copy_matrix_to_tuple_vector!(tuple_vector::AbstractVector{NTuple{N,T}}, matrix::AbstractMatrix{T2}) where {N, T, T2}
    n = Val{N}()
    @parallel for k in axes(matrix,2)
        @inbounds tuple_vector[k] = ntuple(i->(convert(T, @inbounds(matrix[i,k]))), n)
    end
    return tuple_vector
end

for T in (Int32, Int64, Float32, Float64)
    for N in 2:12
        precompile(copy_matrix_to_tuple_vector!, (Vector{NTuple{N, T}}, Matrix{T}))
    end
end

function VoronoiMeshes.PlanarVoronoiDiagram(::Val{N}, nEdges::Vector{UInt8}, ncfile::NCDatasets.NCDataset) where N
    verticesOnCellArray = (ncfile["verticesOnCell"][:,:])::Matrix{Int32}
    nCells = length(nEdges)

    verticesOnCell = ImmutableVectorArray(Vector{NTuple{N, Int32}}(undef,nCells), nEdges)
    t1 = Threads.@spawn copy_matrix_to_tuple_vector!($(verticesOnCell.data), $verticesOnCellArray)

    cellsOnVertexArray = ncfile["cellsOnVertex"][:,:]::Matrix{Int32}
    nVertices = size(cellsOnVertexArray,2)
    cellsOnVertex = Vector{NTuple{3, Int32}}(undef,nVertices)
    t2 = Threads.@spawn copy_matrix_to_tuple_vector!($cellsOnVertex, $cellsOnVertexArray)

    xCell = ncfile["xCell"][:]::Vector{Float64}
    yCell = ncfile["yCell"][:]::Vector{Float64}
    generators = VecArray(x=xCell, y=yCell)

    xVertex = ncfile["xVertex"][:]::Vector{Float64}
    yVertex = ncfile["yVertex"][:]::Vector{Float64}
    vertices = VecArray(x=xVertex, y=yVertex)

    x_period = ncfile.attrib["x_period"]::Float64
    y_period = ncfile.attrib["y_period"]::Float64

    meshDensity = ncfile["meshDensity"][:]::Vector{Float64}

    wait(t1)
    wait(t2)

    return PlanarVoronoiDiagram(generators, vertices, verticesOnCell, cellsOnVertex, meshDensity, x_period, y_period)
end

for N in 6:10
    precompile(VoronoiMeshes.PlanarVoronoiDiagram, (Val{N}, Vector{UInt8}, NCDatasets.NCDataset{Nothing, Missing}))
end

function VoronoiMeshes.PlanarVoronoiDiagram(ncfile::NCDatasets.NCDataset)
    on_a_sphere(ncfile) && throw(error("Mesh is not planar"))

    nEdgesOnCell = UInt8.(ncfile["nEdgesOnCell"][:]::Vector{Int32})
    maxEdges = Int(maximum(nEdgesOnCell))

    #Avoid dynamic dispatch for most common cases
    if maxEdges == 6
        return PlanarVoronoiDiagram(Val{6}(), nEdgesOnCell,ncfile)
    elseif maxEdges == 7
        return PlanarVoronoiDiagram(Val{7}(), nEdgesOnCell,ncfile)
    elseif maxEdges == 8
        return PlanarVoronoiDiagram(Val{8}(), nEdgesOnCell,ncfile)
    elseif maxEdges == 9
        return PlanarVoronoiDiagram(Val{9}(), nEdgesOnCell,ncfile)
    elseif maxEdges == 10
        return PlanarVoronoiDiagram(Val{10}(), nEdgesOnCell,ncfile)
    else
        return PlanarVoronoiDiagram(Val(maxEdges), nEdgesOnCell, ncfile)
    end

end

function VoronoiMeshes.PlanarVoronoiDiagram(n::Val{maxEdges}, ncfile::NCDatasets.NCDataset) where {maxEdges}
    on_a_sphere(ncfile) && throw(error("Mesh is not planar"))

    nEdgesOnCell = UInt8.(ncfile["nEdgesOnCell"][:]::Vector{Int32})
    maxEdges == Int(maximum(nEdgesOnCell)) || throw(error("nEdges not consistent with data in NetCDF file"))

       return PlanarVoronoiDiagram(n, nEdgesOnCell, ncfile)
end

for N in 6:10
    precompile(VoronoiMeshes.PlanarVoronoiDiagram, (Val{N}, NCDatasets.NCDataset{Nothing, Missing}))
end

function VoronoiMeshes.SphericalVoronoiDiagram(::Val{N}, nEdges::Vector{UInt8}, ncfile::NCDatasets.NCDataset) where {N}
    verticesOnCellArray = (ncfile["verticesOnCell"][:,:])::Matrix{Int32}
    nCells = length(nEdges)

    verticesOnCell = ImmutableVectorArray(Vector{NTuple{N, Int32}}(undef,nCells),nEdges)
    t1 = Threads.@spawn copy_matrix_to_tuple_vector!($(verticesOnCell.data), $verticesOnCellArray)

    cellsOnVertexArray = ncfile["cellsOnVertex"][:,:]::Matrix{Int32}
    nVertices = size(cellsOnVertexArray,2)
    cellsOnVertex = Vector{NTuple{3, Int32}}(undef,nVertices)
    t2 = Threads.@spawn copy_matrix_to_tuple_vector!($cellsOnVertex, $cellsOnVertexArray)

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
    precompile(VoronoiMeshes.SphericalVoronoiDiagram, (Val{N}, Vector{UInt8}, NCDatasets.NCDataset{Nothing, Missing}))
end

function VoronoiMeshes.SphericalVoronoiDiagram(ncfile::NCDatasets.NCDataset)
    on_a_sphere(ncfile) || throw(error("Mesh is not spherical"))

    nEdgesOnCell = UInt8.(ncfile["nEdgesOnCell"][:]::Vector{Int32})
    maxEdges = Int(maximum(nEdgesOnCell))

    #Avoid dynamic dispatch for most common cases
    if maxEdges == 6
        return SphericalVoronoiDiagram(Val{6}(), nEdgesOnCell,ncfile)
    elseif maxEdges == 7
        return SphericalVoronoiDiagram(Val{7}(), nEdgesOnCell,ncfile)
    elseif maxEdges == 8
        return SphericalVoronoiDiagram(Val{8}(), nEdgesOnCell,ncfile)
    elseif maxEdges == 9
        return SphericalVoronoiDiagram(Val{9}(), nEdgesOnCell,ncfile)
    elseif maxEdges == 10
        return SphericalVoronoiDiagram(Val{10}(), nEdgesOnCell,ncfile)
    else
        return SphericalVoronoiDiagram(Val(maxEdges), nEdgesOnCell, ncfile)
    end

end

function VoronoiMeshes.SphericalVoronoiDiagram(n::Val{maxEdges}, ncfile::NCDatasets.NCDataset) where {maxEdges}
    on_a_sphere(ncfile) || throw(error("Mesh is not spherical"))

    nEdgesOnCell = UInt8.(ncfile["nEdgesOnCell"][:]::Vector{Int32})
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
    nEdges = UInt8.(ncfile["nEdgesOnCell"][:]::Vector{Int32})
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

    if haskey(ncfile,"localVerticalUnitVectors")
        lvuva = ncfile["localVerticalUnitVectors"]
        n_z = lvuva[3,:]
        if S
            n_x  = lvuva[1,:]
            n_y = lvuva[2,:]
            cell_info.verticalUnitVectors = VecArray(n_x, n_y, n_z) 
        else
            cell_info.verticalUnitVectors = VecArray(z = n_z) 
        end
    end

    if haskey(ncfile,"cellTangentPlane")
        ctp = ncfile["cellTangentPlane"]
        tux = ctp[1,1,:]
        tuy = ctp[2,1,:]
        tvx = ctp[1,2,:]
        tvy = ctp[2,2,:]
        if S
            tuz = ctp[3,1,:]
            tvz = ctp[3,2,:]
            cell_info.tangentPlane = (VecArray(tux, tuy, tuz), VecArray(tvx, tvy, tvz))
        else
            cell_info.tangentPlane = (VecArray(x = tux, y =tuy), VecArray(x = tvx, y = tvy))
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
    edges = ImmutableVectorArray(similar(vertices.data), nEdges)
    copy_matrix_to_tuple_vector!(edges.data, ncfile["edgesOnCell"][:,:]::Matrix{Int32})
    cells = ImmutableVectorArray(similar(vertices.data), nEdges)
    copy_matrix_to_tuple_vector!(cells.data, ncfile["cellsOnCell"][:,:]::Matrix{Int32})

    return Cells(n, position, nEdges, vertices, edges, cells, CellInfo(voro, ncfile))
end

VoronoiMeshes.Cells(ncfile::NCDatasets.NCDataset) = Cells(VoronoiDiagram(ncfile), ncfile)
VoronoiMeshes.Cells(v::Val{N}, ncfile::NCDatasets.NCDataset) where {N} = Cells(VoronoiDiagram(v, ncfile), ncfile)

function VoronoiMeshes.VertexInfo(voro::VoronoiDiagram{S, NE, TI, TF, Tz}, ncfile::NCDatasets.NCDataset) where {S, NE, TI, TF, Tz}
  vertex_info = VertexInfo(voro)

    haskey(ncfile, "centroidVertex") && (vertex_info.area = ncfile["centroidVertex"][:]::Vector{Float64})
    haskey(ncfile, "centroidTriangle") && (vertex_info.area = ncfile["centroidTriangle"][:]::Vector{Float64})
    haskey(ncfile, "areaVertex") && (vertex_info.area = ncfile["areaVertex"][:]::Vector{Float64})
    haskey(ncfile, "areaTriangle") && (vertex_info.area = ncfile["areaTriangle"][:]::Vector{Float64})
    haskey(ncfile, "lonVertex") && (vertex_info.longitude = ncfile["lonVertex"][:]::Vector{Float64})
    haskey(ncfile, "latVertex") && (vertex_info.latitude = ncfile["latVertex"][:]::Vector{Float64})

    if haskey(ncfile, "kiteAreasOnVertex")
        kiteAreasArray = ncfile["kiteAreasOnVertex"][:,:]
        kiteAreas = [(kiteAreasArray[1,k],kiteAreasArray[2,k],kiteAreasArray[3,k]) for k in axes(kiteAreasArray,2)]
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
    copy_matrix_to_tuple_vector!(edges, ncfile["edgesOnVertex"][:,:]::Matrix{Int32})
    return Vertices(n, position, edges, cells, VertexInfo(voro, ncfile))
end

VoronoiMeshes.Vertices(ncfile::NCDatasets.NCDataset) = Vertices(VoronoiDiagram(ncfile), ncfile)
VoronoiMeshes.Vertices(v::Val{N}, ncfile::NCDatasets.NCDataset) where {N} = Vertices(VoronoiDiagram(v, ncfile), ncfile)

for S in (true, false)
    for NE in 6:12
        for TI in (Int32, Int64)
            for TF in (Float32, Float64)
                for Tz in (TF, Zero)
                    precompile(VoronoiMeshes.CellInfo, (VoronoiDiagram{S, NE, TI, TF, Tz}, NCDatasets.NCDataset{Nothing, Missing}))
                    precompile(VoronoiMeshes.Cells, (VoronoiDiagram{S, NE, TI, TF, Tz}, NCDatasets.NCDataset{Nothing, Missing}))
                    precompile(VoronoiMeshes.VertexInfo, (VoronoiDiagram{S, NE, TI, TF, Tz}, NCDatasets.NCDataset{Nothing, Missing}))
                    precompile(VoronoiMeshes.Vertices, (VoronoiDiagram{S, NE, TI, TF, Tz}, NCDatasets.NCDataset{Nothing, Missing}))
                end
            end
        end
    end
end

for func in  (:PlanarVoronoiDiagram, :SphericalVoronoiDiagram, :VoronoiDiagram,
             :CellInfo, :Cells, :VertexInfo, :Vertices)
    @eval begin
        VoronoiMeshes.$func(file_name::String) = NCDataset(file_name) do f; VoronoiMeshes.$func(f);end
        VoronoiMeshes.$func(v::Val, file_name::String) = NCDataset(file_name) do f; VoronoiMeshes.$func(v, f);end
        precompile(VoronoiMeshes.$func,(NCDatasets.NCDataset{Nothing, Missing},))
        precompile(VoronoiMeshes.$func,(String,))
    end
    for N in 6:10
        @eval precompile(VoronoiMeshes.$func,(Val{$N}, NCDatasets.NCDataset{Nothing, Missing},))
    end
end

end # module
#function construct_elementsOnElement!(elementsOnElement::AbstractVector{ImmutableVector{maxEdges,TE}}, elementsOnElementArray::AbstractMatrix{TE}, nElemtensOnElement::AbstractVector{TI}) where {maxEdges,TE,TI}
#    n = Val{maxEdges}()
#    @parallel for k in axes(elementsOnElementArray,2)
#        @inbounds elementsOnElement[k] = ImmutableVector{maxEdges}(ntuple(i->(@inbounds elementsOnElementArray[i,k]), n),nElemtensOnElement[k])
#    end
#    return elementsOnElement
#end
#
#
#function construct_elementsOnElement(n::Val{maxEdges},elementsOnElementArray::AbstractMatrix{TE},nElemtensOnElement::AbstractVector{TI}) where {maxEdges,TE,TI}
#    elementsOnElement = Vector{ImmutableVector{maxEdges,TE}}(undef,size(elementsOnElementArray,2))
#    return construct_elementsOnElement!(elementsOnElement, elementsOnElementArray, nElemtensOnElement)
#end
#
#function VoronoiMeshDataStruct.CellConnectivity(me::Val{maxEdges},nEdgesOnCell::AbstractVector{TI},ncfile::NCDatasets.NCDataset) where {maxEdges,TI}
#    verticesOnCellArray = (ncfile["verticesOnCell"][:,:])::Matrix{Int32}
#    l = UInt8.(nEdgesOnCell)
#    nCells = length(l)
#    verticesOnCell = ImmutableVectorArray(Vector{NTuple{maxEdges, TI}}(undef,nCells),l)
#    t1 = Threads.@spawn copy_matrix_to_tuple_vector!($(verticesOnCell.data), $verticesOnCellArray)
#
#    edgesOnCellArray = (ncfile["edgesOnCell"][:,:])::Matrix{Int32}
#    edgesOnCell = ImmutableVectorArray(Vector{NTuple{maxEdges, TI}}(undef,nCells),l)
#    t2 = Threads.@spawn copy_matrix_to_tuple_vector!($(edgesOnCell.data), $edgesOnCellArray)
#
#    cellsOnCellArray = (ncfile["cellsOnCell"][:,:])::Matrix{Int32}
#    cellsOnCell = ImmutableVectorArray(Vector{NTuple{maxEdges, TI}}(undef,nCells),l)
#    copy_matrix_to_tuple_vector!(cellsOnCell.data, cellsOnCellArray)
#
#    wait(t1)
#    wait(t2)
#    return CellConnectivity(verticesOnCell, cellsOnCell, edgesOnCell)
#end
#
#for N in 6:12
#    for T in (Int64,Int32)
#        for TE in (Int64,Int32,Float64,Float32)
#            precompile(construct_elementsOnElement,(Val{N},Matrix{TE},Vector{T}))
#            precompile(copy_matrix_to_tuple_vector!,(ImVecArray{N, T, 1}, Matrix{TE}))
#        end
#        precompile(VoronoiMeshDataStruct.CellConnectivity,(Val{N},Vector{T},NCDatasets.NCDataset{Nothing,Missing}))
#    end
#end
#
#function VoronoiMeshDataStruct.CellConnectivity(ncfile::NCDatasets.NCDataset)
#    nEdgesOnCell = ncfile["nEdgesOnCell"][:]::Vector{Int32}
#    maxEdges = Int(maximum(nEdgesOnCell))
#
#    #Avoid dynamic dispatch for most common cases
#    if maxEdges == 6
#        return CellConnectivity(Val{6}(), nEdgesOnCell,ncfile)
#    elseif maxEdges == 7
#        return CellConnectivity(Val{7}(), nEdgesOnCell,ncfile)
#    elseif maxEdges == 8
#        return CellConnectivity(Val{8}(), nEdgesOnCell,ncfile)
#    elseif maxEdges == 9
#        return CellConnectivity(Val{9}(), nEdgesOnCell,ncfile)
#    elseif maxEdges == 10
#        return CellConnectivity(Val{10}(), nEdgesOnCell,ncfile)
#    else
#        return CellConnectivity(Val(maxEdges), nEdgesOnCell, ncfile)
#    end
#end
#
#precompile(VoronoiMeshDataStruct.CellConnectivity,(NCDatasets.NCDataset{Nothing,Missing},))
#
#function VoronoiMeshDataStruct.CellBase(onSphere::Val{on_sphere},mE::Val{maxEdges},nEdges,ncfile::NCDatasets.NCDataset) where {on_sphere,maxEdges}
#    indices = CellConnectivity(mE,nEdges,ncfile)
#    x = (ncfile["xCell"][:])::Vector{Float64}
#    y = (ncfile["yCell"][:])::Vector{Float64}
#
#    if on_sphere
#        position3D = VecArray(x = x, y = y, z=(ncfile["zCell"][:])::Vector{Float64})
#        sphere_radius = ncfile.attrib["sphere_radius"]::Float64
#        return CellBase(position3D, indices, sphere_radius)
#    else
#        position2D = VecArray(x = x, y = y)
#        xp = ncfile.attrib["x_period"]::Float64
#        yp = ncfile.attrib["y_period"]::Float64
#        return CellBase(position2D, indices, xp, yp)
#    end
#end
#
#for N in 6:12
#    for on_sphere in (true,false)
#        precompile(VoronoiMeshDataStruct.CellBase,(Val{on_sphere}, Val{N}, Vector{Int32}, NCDatasets.NCDataset{Nothing,Missing}))
#    end
#end
#
#function VoronoiMeshDataStruct.CellBase(ncfile::NCDatasets.NCDataset)
#    nEdges = ncfile["nEdgesOnCell"][:]::Vector{Int32}
#    maxEdges = Int(maximum(nEdges))
#    onSphere, on_sphere = _on_a_sphere(ncfile)
#
#    #Avoid dynamic dispatch for most common cases
#    if on_sphere
#        if maxEdges == 6
#            return CellBase(Val{true}(), Val{6}(), nEdges, ncfile)
#        elseif maxEdges == 7
#            return CellBase(Val{true}(), Val{7}(), nEdges, ncfile)
#        elseif maxEdges == 8
#            return CellBase(Val{true}(), Val{8}(), nEdges, ncfile)
#        elseif maxEdges == 9
#            return CellBase(Val{true}(), Val{9}(), nEdges, ncfile)
#        elseif maxEdges == 10
#            return CellBase(Val{true}(), Val{10}(), nEdges, ncfile)
#        end
#    else
#        if maxEdges == 6
#            return CellBase(Val{false}(), Val{6}(), nEdges, ncfile)
#        elseif maxEdges == 7
#            return CellBase(Val{false}(), Val{7}(), nEdges, ncfile)
#        elseif maxEdges == 8
#            return CellBase(Val{false}(), Val{8}(), nEdges, ncfile)
#        elseif maxEdges == 9
#            return CellBase(Val{false}(), Val{9}(), nEdges, ncfile)
#        elseif maxEdges == 10
#            return CellBase(Val{false}(), Val{10}(), nEdges, ncfile)
#        end
#    end
#    return CellBase(onSphere,Val(maxEdges),nEdges,ncfile)
#end
#
#precompile(VoronoiMeshDataStruct.CellBase,(NCDatasets.NCDataset{Nothing,Missing},))
#
##const cell_info_vectors = (longitude="lonCell", latitude="latCell",
##                          meshDensity="meshDensity",indexToID="indexToCellID",
##                          area="areaCell", bdyMask="bdyMaskCell")
#
#const cell_info_vectors = (longitude="lonCell", latitude="latCell",
#                          meshDensity="meshDensity", area="areaCell")
#
##const cell_info_matrices_max_edges = (defcA="defc_a", defcB="defc_b",
##                                      xGradientCoeff="cell_gradient_coef_x", yGradientCoeff="cell_gradient_coef_y")
#
#function VoronoiMeshDataStruct.Cells(ncfile::NCDatasets.NCDataset)
#    cells = CellBase(ncfile)
#    return Cells(cells,ncfile)
#end
#precompile(VoronoiMeshDataStruct.Cells,(NCDatasets.NCDataset{Nothing,Missing},))
#
#function VoronoiMeshDataStruct.Cells(cells::CellBase{on_sphere,max_nedges},ncfile::NCDatasets.NCDataset) where {on_sphere,max_nedges}
#
#    longitude = ncfile["lonCell"][:]::Vector{Float64}
#    latitude = ncfile["latCell"][:]::Vector{Float64}
#    area = ncfile["areaCell"][:]::Vector{Float64}
#    meshDensity = ncfile["meshDensity"][:]::Vector{Float64}
#
#    if haskey(ncfile,"localVerticalUnitVectors")
#        lvuva = ncfile["localVerticalUnitVectors"]
#        verticalUnitVectors = on_sphere ? VecArray(x=lvuva[1,:], y=lvuva[2,:], z=lvuva[3,:]) : VecArray(z=lvuva[3,:])
#    else
#        verticalUnitVectors = on_sphere ? normalize(cells.position) : VecArray(z=ones(cells.n))
#    end
#
#    if haskey(ncfile,"cellTangentPlane")
#        ctp = ncfile["cellTangentPlane"]
#        tangentPlane = on_sphere ?
#                                (VecArray(x=ctp[1,1,:],y=ctp[2,1,:],z=ctp[3,1,:]),VecArray(x=ctp[1,2,:],y=ctp[2,2,:],z=ctp[3,2,:])) :
#                                (VecArray(x=ctp[1,1,:],y=ctp[2,1,:]),VecArray(x=ctp[1,2,:],y=ctp[2,2,:]))
#    else
#        tangentPlane = if on_sphere
#            (VecArray(x=zeros(0), y = zeros(0), z= zeros(0)),
#             VecArray(x=zeros(0), y = zeros(0), z= zeros(0)))
#        else
#            (VecArray(x=zeros(0), y = zeros(0)),
#             VecArray(x=zeros(0), y = zeros(0)))
#        end
#    end
#
#    return Cells(cells, longitude, latitude, area, verticalUnitVectors, tangentPlane, meshDensity)
#end
#
#for N in 6:12
#    for on_sphere in (true,false)
#        for TF in (Float64,Float32)
#            for Tz in (Float64,Float32,Zero)
#                precompile(VoronoiMeshDataStruct.Cells,(VoronoiMeshDataStruct.CellBase{on_sphere,N,Int32,TF,Tz}, NCDatasets.NCDataset{Nothing,Missing}))
#            end
#        end
#    end
#end
#
#function VoronoiMeshDataStruct.VertexConnectivity(ncfile::NCDatasets.NCDataset)
#    edgesOnVertexArray = ncfile["edgesOnVertex"][:,:]::Matrix{Int32}
#    edgesOnVertex = [(edgesOnVertexArray[1,k],edgesOnVertexArray[2,k],edgesOnVertexArray[3,k]) for k in axes(edgesOnVertexArray,2)]
#
#    cellsOnVertexArray = ncfile["cellsOnVertex"][:,:]::Matrix{Int32}
#    cellsOnVertex = [(cellsOnVertexArray[1,k],cellsOnVertexArray[2,k],cellsOnVertexArray[3,k]) for k in axes(cellsOnVertexArray,2)]
#    return VertexConnectivity(edgesOnVertex,cellsOnVertex)
#end
#
#precompile(VoronoiMeshDataStruct.VertexConnectivity,(NCDatasets.NCDataset{Nothing,Missing},))
#
#function VoronoiMeshDataStruct.VertexBase(ncfile::NCDatasets.NCDataset)
#    indices = VertexConnectivity(ncfile)
#
#    _, on_sphere = _on_a_sphere(ncfile)
#    
#    x = (ncfile["xVertex"][:])::Vector{Float64}
#    y = (ncfile["yVertex"][:])::Vector{Float64}
#
#    if on_sphere
#        position = VecArray(x = x, y = y, z = (ncfile["zVertex"][:])::Vector{Float64})
#        return VertexBase(length(position.x), indices, position, Val{true}())
#    else
#        position_p = VecArray(x = x, y = y)
#        return VertexBase(length(position_p.x), indices, position_p, Val{false}())
#    end
#end
#
#precompile(VoronoiMeshDataStruct.VertexBase,(NCDatasets.NCDataset{Nothing,Missing},))
#
#const vertex_info_vectors = (longitude="lonVertex", latitude="latVertex",
#                             indexToID="indexToVertexID",
#                             area="areaTriangle", bdyMask="bdyMaskVertex")
#
#function VoronoiMeshDataStruct.VertexInfo(ncfile::NCDatasets.NCDataset)
#    vertexBase = VertexBase(ncfile)
#    vertex = VertexInfo(vertexBase)
#
#    for (field_name, nc_name) in pairs(vertex_info_vectors)
#        if haskey(ncfile,nc_name)
#            setproperty!(vertex,field_name,ncfile[nc_name][:])
#        end
#    end
#
#    if haskey(ncfile,"kiteAreasOnVertex")
#        kiteAreasArray = ncfile["kiteAreasOnVertex"]
#        kiteAreas = [(kiteAreasArray[1,k],kiteAreasArray[2,k],kiteAreasArray[3,k]) for k in axes(kiteAreasArray,2)]
#        vertex.kiteAreas = kiteAreas
#    end
#
#    return vertex
#end
#
#precompile(VoronoiMeshDataStruct.VertexInfo,(NCDatasets.NCDataset{Nothing,Missing},))
#
#function VoronoiMeshDataStruct.EdgeConnectivity(ncfile::NCDatasets.NCDataset)
#    verticesOnEdgeArray = ncfile["verticesOnEdge"][:,:]::Matrix{Int32}
#    verticesOnEdge = [(verticesOnEdgeArray[1,k],verticesOnEdgeArray[2,k]) for k in axes(verticesOnEdgeArray,2)]
#
#    cellsOnEdgeArray = ncfile["cellsOnEdge"][:,:]::Matrix{Int32}
#    cellsOnEdge = [(cellsOnEdgeArray[1,k],cellsOnEdgeArray[2,k]) for k in axes(cellsOnEdgeArray,2)]
#    return EdgeConnectivity(verticesOnEdge,cellsOnEdge)
#end
#
#precompile(VoronoiMeshDataStruct.EdgeConnectivity,(NCDatasets.NCDataset{Nothing,Missing},))
#
#function VoronoiMeshDataStruct.EdgeBase(ncfile::NCDatasets.NCDataset)
#    indices = EdgeConnectivity(ncfile)
#
#    _, on_sphere = _on_a_sphere(ncfile)
#
#    x = (ncfile["xEdge"][:])::Vector{Float64}
#    y = (ncfile["yEdge"][:])::Vector{Float64}
#
#    if on_sphere
#        position = VecArray(x = x, y = y, z = (ncfile["zEdge"][:])::Vector{Float64})
#        return EdgeBase(length(x), indices, position, Val{true}())
#    else
#        position_p = VecArray(x = x, y = y)
#        return EdgeBase(length(x), indices, position_p, Val{false}())
#    end
#end
#
#precompile(VoronoiMeshDataStruct.EdgeBase,(NCDatasets.NCDataset{Nothing,Missing},))
#
#function VoronoiMeshDataStruct.EdgeVelocityReconstruction(ncfile::NCDatasets.NCDataset)
#    nEdges = ncfile["nEdgesOnEdge"][:]::Vector{Int32}
#    max_n_edges = Int(maximum(nEdges))
#    return EdgeVelocityReconstruction(Val(max_n_edges),nEdges,ncfile)
#end
#
#precompile(VoronoiMeshDataStruct.EdgeVelocityReconstruction,(NCDatasets.NCDataset{Nothing,Missing},))
#
#function VoronoiMeshDataStruct.EdgeVelocityReconstruction(ne::Val{max_n_edges},nEdges,ncfile::NCDatasets.NCDataset) where {max_n_edges}
#    edgesOnEdgeArray = ncfile["edgesOnEdge"][:,:]::Matrix{Int32}
#    indices = construct_elementsOnElement(ne,edgesOnEdgeArray,nEdges)
#    weightsArray = ncfile["weightsOnEdge"][:,:]
#    weights = construct_elementsOnElement(ne,weightsArray,nEdges)
#    return EdgeVelocityReconstruction(nEdges,indices,weights)
#end
#
#for N in 10:24
#    precompile(VoronoiMeshDataStruct.EdgeVelocityReconstruction,(Val{N},Vector{Int32},NCDatasets.NCDataset{Nothing,Missing}))
#end
#
#const edge_info_vectors = (longitude="lonEdge", latitude="latEdge",
#                           indexToID="indexToEdgeID", dv="dvEdge", dc="dcEdge",
#                           angle="angleEdge", bdyMask="bdyMaskEdge")
#
#function VoronoiMeshDataStruct.EdgeInfo(ncfile::NCDatasets.NCDataset)
#    edgeinfo = EdgeInfo(EdgeBase(ncfile),EdgeVelocityReconstruction(ncfile))
#
#    _, on_sphere = _on_a_sphere(ncfile)
#
#    for (field_name, nc_name) in pairs(edge_info_vectors)
#        if haskey(ncfile,nc_name)
#            setproperty!(edgeinfo,field_name,ncfile[nc_name][:])
#        end
#    end
#
#    if haskey(ncfile,"edgeNormalVectors")
#        env = ncfile["edgeNormalVectors"]
#        edgeinfo.normalVectors = on_sphere ? VecArray(x=env[1,:], y=env[2,:], z=env[3,:]) : VecArray(x=env[1,:], y=env[2,:])
#    end
#
#    if haskey(ncfile,"deriv_two")
#        edgeinfo.derivTwo = ncfile["deriv_two"][:,:,:]
#    end
#
#    return edgeinfo
#end
#
#precompile(VoronoiMeshDataStruct.EdgeInfo,(NCDatasets.NCDataset{Nothing,Missing},))
#
#function VoronoiMeshDataStruct.VoronoiMesh(ncfile::NCDatasets.NCDataset)
#    attributes = Dict{Symbol,Union{String,Float64,Float32,Int64,Int32}}()
#    for (key,val) in ncfile.attrib
#        val isa String && (val = String(strip(val)))
#        attributes[Symbol(key)] = val
#    end
#    return VoronoiMesh(Cells(ncfile), VertexInfo(ncfile), EdgeInfo(ncfile), attributes) 
#end
#
#precompile(VoronoiMeshDataStruct.VoronoiMesh,(NCDatasets.NCDataset{Nothing,Missing},))
#
#for func in  (:VoronoiMesh,
#              :Cells,:CellBase,:CellConnectivity,
#              :VertexInfo,:VertexBase,:VertexConnectivity,
#              :EdgeInfo,:EdgeBase,:EdgeConnectivity,:EdgeVelocityReconstruction)
#
#    @eval begin
#        VoronoiMeshDataStruct.$func(file_name::String) = NCDataset(file_name) do f; VoronoiMeshDataStruct.$func(f);end
#        precompile(VoronoiMeshDataStruct.$func,(String,))
#    end
#end

