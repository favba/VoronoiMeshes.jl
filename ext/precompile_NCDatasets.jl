precompile(Tuple{typeof(NCDatasets.nc_open), String, UInt16})
precompile(Tuple{Type{NCDatasets.Variable{Int32, 1, NCDatasets.NCDataset{Nothing, Base.Missing}}}, NCDatasets.NCDataset{Nothing, Base.Missing}, Int32, Tuple{Int32}})
precompile(Tuple{Type{NCDatasets.CommonDataModel.Attributes{TDS} where TDS<:Union{NCDatasets.CommonDataModel.AbstractDataset, NCDatasets.CommonDataModel.AbstractVariable{T, N} where N where T}}, NCDatasets.Variable{Int32, 1, NCDatasets.NCDataset{Nothing, Base.Missing}}})
precompile(Tuple{Type{NCDatasets.Variable{Int32, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}}, NCDatasets.NCDataset{Nothing, Base.Missing}, Int32, Tuple{Int32, Int32}})
precompile(Tuple{Type{NCDatasets.CommonDataModel.Attributes{TDS} where TDS<:Union{NCDatasets.CommonDataModel.AbstractDataset, NCDatasets.CommonDataModel.AbstractVariable{T, N} where N where T}}, NCDatasets.Variable{Int32, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}})
precompile(Tuple{Type{NCDatasets.Variable{Float64, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}}, NCDatasets.NCDataset{Nothing, Base.Missing}, Int32, Tuple{Int32, Int32}})
precompile(Tuple{Type{NCDatasets.CommonDataModel.Attributes{TDS} where TDS<:Union{NCDatasets.CommonDataModel.AbstractDataset, NCDatasets.CommonDataModel.AbstractVariable{T, N} where N where T}}, NCDatasets.Variable{Float64, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}})
precompile(Tuple{Type{NCDatasets.Variable{Float64, 1, NCDatasets.NCDataset{Nothing, Base.Missing}}}, NCDatasets.NCDataset{Nothing, Base.Missing}, Int32, Tuple{Int32}})
precompile(Tuple{Type{NCDatasets.CommonDataModel.Attributes{TDS} where TDS<:Union{NCDatasets.CommonDataModel.AbstractDataset, NCDatasets.CommonDataModel.AbstractVariable{T, N} where N where T}}, NCDatasets.Variable{Float64, 1, NCDatasets.NCDataset{Nothing, Base.Missing}}})
precompile(Tuple{typeof(NCDatasets.CommonDataModel._getattrib), NCDatasets.NCDataset{Nothing, Base.Missing}, NCDatasets.Variable{Int32, 1, NCDatasets.NCDataset{Nothing, Base.Missing}}, String, String, Nothing})
precompile(Tuple{Type{NCDatasets.CommonDataModel.CFVariable{Int32, 1, NCDatasets.Variable{Int32, 1, NCDatasets.NCDataset{Nothing, Base.Missing}}, NCDatasets.CommonDataModel.Attributes{NCDatasets.Variable{Int32, 1, NCDatasets.NCDataset{Nothing, Base.Missing}}}, NamedTuple{(:fillvalue, :missing_values, :scale_factor, :add_offset, :calendar, :time_origin, :time_factor, :maskingvalue), Tuple{Nothing, Tuple{}, Nothing, Nothing, Nothing, Nothing, Nothing, Base.Missing}}}}, NCDatasets.Variable{Int32, 1, NCDatasets.NCDataset{Nothing, Base.Missing}}, NCDatasets.CommonDataModel.Attributes{NCDatasets.Variable{Int32, 1, NCDatasets.NCDataset{Nothing, Base.Missing}}}, NamedTuple{(:fillvalue, :missing_values, :scale_factor, :add_offset, :calendar, :time_origin, :time_factor, :maskingvalue), Tuple{Nothing, Tuple{}, Nothing, Nothing, Nothing, Nothing, Nothing, Base.Missing}}})
precompile(Tuple{typeof(Base.getproperty), NCDatasets.DiskArrays.GridChunks{0, Tuple{}}, Symbol})
precompile(Tuple{Type{NCDatasets.DiskArrays.SubRanges{S} where S}, NCDatasets.DiskArrays.CanStepRange, Float64})
precompile(Tuple{Type{NCDatasets.DiskArrays.Unchunked{BS} where BS}, NCDatasets.DiskArrays.SubRanges{NCDatasets.DiskArrays.CanStepRange}})
precompile(Tuple{Type{NCDatasets.DiskArrays.Chunked{BS} where BS}, NCDatasets.DiskArrays.SubRanges{NCDatasets.DiskArrays.CanStepRange}})
precompile(Tuple{Type{NCDatasets.DiskArrays.NoBatch{S} where S}, Bool, Float64})
precompile(Tuple{typeof(Base.getproperty), NamedTuple{(:fillvalue, :missing_values, :scale_factor, :add_offset, :calendar, :time_origin, :time_factor, :maskingvalue), Tuple{Nothing, Tuple{}, Nothing, Nothing, Nothing, Nothing, Nothing, Base.Missing}}, Symbol})
precompile(Tuple{typeof(Base.getindex), NCDatasets.CommonDataModel.CFVariable{Int32, 1, NCDatasets.Variable{Int32, 1, NCDatasets.NCDataset{Nothing, Base.Missing}}, NCDatasets.CommonDataModel.Attributes{NCDatasets.Variable{Int32, 1, NCDatasets.NCDataset{Nothing, Base.Missing}}}, NamedTuple{(:fillvalue, :missing_values, :scale_factor, :add_offset, :calendar, :time_origin, :time_factor, :maskingvalue), Tuple{Nothing, Tuple{}, Nothing, Nothing, Nothing, Nothing, Nothing, Base.Missing}}}, Base.Colon})
precompile(Tuple{typeof(Base.reverse), Tuple{Int64}})
precompile(Tuple{typeof(NCDatasets.CommonDataModel._getattrib), NCDatasets.NCDataset{Nothing, Base.Missing}, NCDatasets.Variable{Int32, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, String, String, Nothing})
precompile(Tuple{Type{NCDatasets.CommonDataModel.CFVariable{Int32, 2, NCDatasets.Variable{Int32, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, NCDatasets.CommonDataModel.Attributes{NCDatasets.Variable{Int32, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}}, NamedTuple{(:fillvalue, :missing_values, :scale_factor, :add_offset, :calendar, :time_origin, :time_factor, :maskingvalue), Tuple{Nothing, Tuple{}, Nothing, Nothing, Nothing, Nothing, Nothing, Base.Missing}}}}, NCDatasets.Variable{Int32, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, NCDatasets.CommonDataModel.Attributes{NCDatasets.Variable{Int32, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}}, NamedTuple{(:fillvalue, :missing_values, :scale_factor, :add_offset, :calendar, :time_origin, :time_factor, :maskingvalue), Tuple{Nothing, Tuple{}, Nothing, Nothing, Nothing, Nothing, Nothing, Base.Missing}}})
precompile(Tuple{typeof(Base.getproperty), NCDatasets.DiskArrays.DiskIndex{0, 0, Tuple{}, Tuple{}, Tuple{}}, Symbol})
precompile(Tuple{typeof(Base.getindex), NCDatasets.CommonDataModel.CFVariable{Int32, 2, NCDatasets.Variable{Int32, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, NCDatasets.CommonDataModel.Attributes{NCDatasets.Variable{Int32, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}}, NamedTuple{(:fillvalue, :missing_values, :scale_factor, :add_offset, :calendar, :time_origin, :time_factor, :maskingvalue), Tuple{Nothing, Tuple{}, Nothing, Nothing, Nothing, Nothing, Nothing, Base.Missing}}}, Base.Colon, Base.Colon})
precompile(Tuple{typeof(Base.reverse), Tuple{Int64, Int64}})
precompile(Tuple{typeof(NCDatasets.size_getindex), NCDatasets.Variable{Int32, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, Base.UnitRange{Int64}, Vararg{Base.UnitRange{Int64}}})
precompile(Tuple{typeof(NCDatasets._size_getindex), NCDatasets.Variable{Int32, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, Tuple{}, Int64, Base.UnitRange{Int64}, Base.UnitRange{Int64}})
precompile(Tuple{typeof(NCDatasets.CommonDataModel._getattrib), NCDatasets.NCDataset{Nothing, Base.Missing}, NCDatasets.Variable{Float64, 1, NCDatasets.NCDataset{Nothing, Base.Missing}}, String, String, Nothing})
precompile(Tuple{Type{NCDatasets.CommonDataModel.CFVariable{Float64, 1, NCDatasets.Variable{Float64, 1, NCDatasets.NCDataset{Nothing, Base.Missing}}, NCDatasets.CommonDataModel.Attributes{NCDatasets.Variable{Float64, 1, NCDatasets.NCDataset{Nothing, Base.Missing}}}, NamedTuple{(:fillvalue, :missing_values, :scale_factor, :add_offset, :calendar, :time_origin, :time_factor, :maskingvalue), Tuple{Nothing, Tuple{}, Nothing, Nothing, Nothing, Nothing, Nothing, Base.Missing}}}}, NCDatasets.Variable{Float64, 1, NCDatasets.NCDataset{Nothing, Base.Missing}}, NCDatasets.CommonDataModel.Attributes{NCDatasets.Variable{Float64, 1, NCDatasets.NCDataset{Nothing, Base.Missing}}}, NamedTuple{(:fillvalue, :missing_values, :scale_factor, :add_offset, :calendar, :time_origin, :time_factor, :maskingvalue), Tuple{Nothing, Tuple{}, Nothing, Nothing, Nothing, Nothing, Nothing, Base.Missing}}})
precompile(Tuple{typeof(Base.getindex), NCDatasets.CommonDataModel.CFVariable{Float64, 1, NCDatasets.Variable{Float64, 1, NCDatasets.NCDataset{Nothing, Base.Missing}}, NCDatasets.CommonDataModel.Attributes{NCDatasets.Variable{Float64, 1, NCDatasets.NCDataset{Nothing, Base.Missing}}}, NamedTuple{(:fillvalue, :missing_values, :scale_factor, :add_offset, :calendar, :time_origin, :time_factor, :maskingvalue), Tuple{Nothing, Tuple{}, Nothing, Nothing, Nothing, Nothing, Nothing, Base.Missing}}}, Base.Colon})
precompile(Tuple{typeof(NCDatasets.CommonDataModel._getattrib), NCDatasets.NCDataset{Nothing, Base.Missing}, NCDatasets.Variable{Float64, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, String, String, Nothing})
precompile(Tuple{Type{NCDatasets.CommonDataModel.CFVariable{Float64, 2, NCDatasets.Variable{Float64, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, NCDatasets.CommonDataModel.Attributes{NCDatasets.Variable{Float64, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}}, NamedTuple{(:fillvalue, :missing_values, :scale_factor, :add_offset, :calendar, :time_origin, :time_factor, :maskingvalue), Tuple{Nothing, Tuple{}, Nothing, Nothing, Nothing, Nothing, Nothing, Base.Missing}}}}, NCDatasets.Variable{Float64, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, NCDatasets.CommonDataModel.Attributes{NCDatasets.Variable{Float64, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}}, NamedTuple{(:fillvalue, :missing_values, :scale_factor, :add_offset, :calendar, :time_origin, :time_factor, :maskingvalue), Tuple{Nothing, Tuple{}, Nothing, Nothing, Nothing, Nothing, Nothing, Base.Missing}}})
precompile(Tuple{typeof(Base.getindex), NCDatasets.CommonDataModel.CFVariable{Float64, 2, NCDatasets.Variable{Float64, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, NCDatasets.CommonDataModel.Attributes{NCDatasets.Variable{Float64, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}}, NamedTuple{(:fillvalue, :missing_values, :scale_factor, :add_offset, :calendar, :time_origin, :time_factor, :maskingvalue), Tuple{Nothing, Tuple{}, Nothing, Nothing, Nothing, Nothing, Nothing, Base.Missing}}}, Base.Colon, Base.Colon})
precompile(Tuple{typeof(NCDatasets.size_getindex), NCDatasets.Variable{Float64, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, Base.UnitRange{Int64}, Vararg{Base.UnitRange{Int64}}})
precompile(Tuple{typeof(NCDatasets._size_getindex), NCDatasets.Variable{Float64, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, Tuple{}, Int64, Base.UnitRange{Int64}, Base.UnitRange{Int64}})
for N in 6:9
    precompile(Tuple{typeof(Base.setproperty!), VoronoiMeshes.VertexInfo{false, N, Int32, Float64, Zeros.Zero}, Symbol, Array{Tuple{Float64, Float64, Float64}, 1}})
    precompile(Tuple{typeof(Base.setproperty!), VoronoiMeshes.VertexInfo{true, N, Int32, Float64, Float64}, Symbol, Array{Tuple{Float64, Float64, Float64}, 1}})
end
precompile(Tuple{Type{NCDatasets.Variable{Float64, 0, NCDatasets.NCDataset{Nothing, Base.Missing}}}, NCDatasets.NCDataset{Nothing, Base.Missing}, Int32, Tuple{}})
precompile(Tuple{Type{NCDatasets.CommonDataModel.Attributes{TDS} where TDS<:Union{NCDatasets.CommonDataModel.AbstractDataset, NCDatasets.CommonDataModel.AbstractVariable{T, N} where N where T}}, NCDatasets.Variable{Float64, 0, NCDatasets.NCDataset{Nothing, Base.Missing}}})
precompile(Tuple{Type{NCDatasets.Variable{Char, 1, NCDatasets.NCDataset{Nothing, Base.Missing}}}, NCDatasets.NCDataset{Nothing, Base.Missing}, Int32, Tuple{Int32}})
precompile(Tuple{Type{NCDatasets.CommonDataModel.Attributes{TDS} where TDS<:Union{NCDatasets.CommonDataModel.AbstractDataset, NCDatasets.CommonDataModel.AbstractVariable{T, N} where N where T}}, NCDatasets.Variable{Char, 1, NCDatasets.NCDataset{Nothing, Base.Missing}}})
precompile(Tuple{Type{NamedTuple{(:force3D, :write_computed), T} where T<:Tuple}, Tuple{Bool, Bool}})
precompile(Tuple{Type{NCDatasets.DiskArrays.DiskIndex{N, M, A, B, C} where C<:Tuple where B<:Tuple where A<:Tuple where M where N}, Tuple{}, Tuple{Int64}, Tuple{}, Tuple{Int64}, Tuple{Base.UnitRange{Int64}}})
precompile(Tuple{typeof(NCDatasets.DiskArrays.merge_index), NCDatasets.DiskArrays.DiskIndex{0, 0, Tuple{}, Tuple{}, Tuple{}}, NCDatasets.DiskArrays.DiskIndex{0, 1, Tuple{}, Tuple{Int64}, Tuple{Base.UnitRange{Int64}}}})
precompile(Tuple{Type{NamedTuple{(:density_threshold,), T} where T<:Tuple}, Tuple{Float64}})
precompile(Tuple{Type{NamedTuple{(:format,), T} where T<:Tuple}, Tuple{Symbol}})
for N in 6:10
    precompile(Tuple{typeof(Core.kwcall), NamedTuple{(:force3D, :write_computed), Tuple{Bool, Bool}}, typeof(VoronoiMeshes.save), String, VoronoiMeshes.VoronoiMesh{false, N, Int32, Float64, Zeros.Zero}})
    precompile(Tuple{typeof(Core.kwcall), NamedTuple{(:force3D, :write_computed), Tuple{Bool, Bool}}, typeof(VoronoiMeshes.save), String, VoronoiMeshes.VoronoiMesh{true, N, Int32, Float64, Float64}})
end
precompile(Tuple{typeof(Core.kwcall), NamedTuple{(:format,), Tuple{Symbol}}, Type{NCDatasets.NCDataset{TDS, Tmaskingvalue} where Tmaskingvalue where TDS}, String, String})
precompile(Tuple{typeof(NCDatasets.nc_create), String, UInt16})
precompile(Tuple{typeof(Base.setindex!), NCDatasets.CommonDataModel.Attributes{NCDatasets.Variable{Float64, 1, NCDatasets.NCDataset{Nothing, Base.Missing}}}, String, String})
@static if pkgversion(NCDatasets.DiskArrays) < v"0.4.10"
    precompile(Tuple{typeof(Base.argtail), NCDatasets.DiskArrays.MRArray{Any, 0, Tuple{}}})
    precompile(Tuple{typeof(Base.getproperty), NCDatasets.DiskArrays.MRArray{Any, 0, Tuple{}}, Symbol})
else
    precompile(Tuple{typeof(Base.argtail), NCDatasets.DiskArrays.MultiReadArray{Any, 0, Tuple{}}})
    precompile(Tuple{typeof(Base.getproperty), NCDatasets.DiskArrays.MultiReadArray{Any, 0, Tuple{}}, Symbol})
end
precompile(Tuple{typeof(Base.setindex!), NCDatasets.CommonDataModel.CFVariable{Float64, 1, NCDatasets.Variable{Float64, 1, NCDatasets.NCDataset{Nothing, Base.Missing}}, NCDatasets.CommonDataModel.Attributes{NCDatasets.Variable{Float64, 1, NCDatasets.NCDataset{Nothing, Base.Missing}}}, NamedTuple{(:fillvalue, :missing_values, :scale_factor, :add_offset, :calendar, :time_origin, :time_factor, :maskingvalue), Tuple{Nothing, Tuple{}, Nothing, Nothing, Nothing, Nothing, Nothing, Base.Missing}}}, Array{Float64, 1}, Base.OneTo{Int64}})
@static if pkgversion(NCDatasets.DiskArrays) < v"0.4.10"
    precompile(Tuple{typeof(NCDatasets.DiskArrays.writeblock_sizecheck!), NCDatasets.Variable{Float64, 1, NCDatasets.NCDataset{Nothing, Base.Missing}}, Array{Float64, 1}, Base.OneTo{Int64}})
else
    precompile(Tuple{typeof(NCDatasets.DiskArrays.writeblock_checked!), NCDatasets.Variable{Float64, 1, NCDatasets.NCDataset{Nothing, Base.Missing}}, Array{Float64, 1}, Base.OneTo{Int64}})
end
precompile(Tuple{typeof(Base.setindex!), NCDatasets.CommonDataModel.Attributes{NCDatasets.Variable{Int32, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}}, String, String})
precompile(Tuple{typeof(Base.setindex!), NCDatasets.CommonDataModel.CFVariable{Int32, 2, NCDatasets.Variable{Int32, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, NCDatasets.CommonDataModel.Attributes{NCDatasets.Variable{Int32, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}}, NamedTuple{(:fillvalue, :missing_values, :scale_factor, :add_offset, :calendar, :time_origin, :time_factor, :maskingvalue), Tuple{Nothing, Tuple{}, Nothing, Nothing, Nothing, Nothing, Nothing, Base.Missing}}}, Base.ReinterpretArray{Int32, 2, Tuple{Int32, Int32, Int32}, Array{Tuple{Int32, Int32, Int32}, 1}, true}, Base.OneTo{Int64}, Base.OneTo{Int64}})
@static if pkgversion(NCDatasets.DiskArrays) < v"0.4.10"
    precompile(Tuple{typeof(NCDatasets.DiskArrays.writeblock_sizecheck!), NCDatasets.Variable{Int32, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, Base.ReinterpretArray{Int32, 2, Tuple{Int32, Int32, Int32}, Array{Tuple{Int32, Int32, Int32}, 1}, true}, Base.OneTo{Int64}, Vararg{Base.OneTo{Int64}}})
else
    precompile(Tuple{typeof(NCDatasets.DiskArrays.writeblock_checked!), NCDatasets.Variable{Int32, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, Base.ReinterpretArray{Int32, 2, Tuple{Int32, Int32, Int32}, Array{Tuple{Int32, Int32, Int32}, 1}, true}, Base.OneTo{Int64}, Vararg{Base.OneTo{Int64}}})
end
precompile(Tuple{typeof(NCDatasets.DiskArrays.writeblock!), NCDatasets.Variable{Int32, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, Base.ReinterpretArray{Int32, 2, Tuple{Int32, Int32, Int32}, Array{Tuple{Int32, Int32, Int32}, 1}, true}, Base.OneTo{Int64}, Vararg{Base.OneTo{Int64}}})
precompile(Tuple{typeof(NCDatasets._write_data_to_nc), NCDatasets.Variable{Int32, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, Base.ReinterpretArray{Int32, 2, Tuple{Int32, Int32, Int32}, Array{Tuple{Int32, Int32, Int32}, 1}, true}, Base.OneTo{Int64}, Vararg{Base.OneTo{Int64}}})
precompile(Tuple{typeof(Base.to_indices), NCDatasets.Variable{Int32, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, Tuple{Base.OneTo{Int64}, Base.OneTo{Int64}}, Tuple{Base.OneTo{Int64}, Base.OneTo{Int64}}})
precompile(Tuple{typeof(NCDatasets.size_getindex), NCDatasets.Variable{Int32, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, Base.OneTo{Int64}, Vararg{Base.OneTo{Int64}}})
precompile(Tuple{typeof(NCDatasets._size_getindex), NCDatasets.Variable{Int32, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, Tuple{}, Int64, Base.OneTo{Int64}, Base.OneTo{Int64}})
for N in 6:9
    precompile(Tuple{typeof(Base.setindex!), NCDatasets.CommonDataModel.CFVariable{Int32, 2, NCDatasets.Variable{Int32, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, NCDatasets.CommonDataModel.Attributes{NCDatasets.Variable{Int32, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}}, NamedTuple{(:fillvalue, :missing_values, :scale_factor, :add_offset, :calendar, :time_origin, :time_factor, :maskingvalue), Tuple{Nothing, Tuple{}, Nothing, Nothing, Nothing, Nothing, Nothing, Base.Missing}}}, Base.ReinterpretArray{Int32, 2, NTuple{N, Int32}, Array{NTuple{N, Int32}, 1}, true}, Base.OneTo{Int64}, Base.OneTo{Int64}})
    @static if pkgversion(NCDatasets.DiskArrays) < v"0.4.10"
        precompile(Tuple{typeof(NCDatasets.DiskArrays.writeblock_sizecheck!), NCDatasets.Variable{Int32, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, Base.ReinterpretArray{Int32, 2, NTuple{N, Int32}, Array{NTuple{N, Int32}, 1}, true}, Base.OneTo{Int64}, Vararg{Base.OneTo{Int64}}})
    else
        precompile(Tuple{typeof(NCDatasets.DiskArrays.writeblock_checked!), NCDatasets.Variable{Int32, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, Base.ReinterpretArray{Int32, 2, NTuple{N, Int32}, Array{NTuple{N, Int32}, 1}, true}, Base.OneTo{Int64}, Vararg{Base.OneTo{Int64}}})
    end
    precompile(Tuple{typeof(NCDatasets.DiskArrays.writeblock!), NCDatasets.Variable{Int32, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, Base.ReinterpretArray{Int32, 2, NTuple{N, Int32}, Array{NTuple{N, Int32}, 1}, true}, Base.OneTo{Int64}, Vararg{Base.OneTo{Int64}}})
    precompile(Tuple{typeof(NCDatasets._write_data_to_nc), NCDatasets.Variable{Int32, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, Base.ReinterpretArray{Int32, 2, NTuple{N, Int32}, Array{NTuple{N, Int32}, 1}, true}, Base.OneTo{Int64}, Vararg{Base.OneTo{Int64}}})
end
precompile(Tuple{typeof(Base.setindex!), NCDatasets.CommonDataModel.Attributes{NCDatasets.Variable{Int32, 1, NCDatasets.NCDataset{Nothing, Base.Missing}}}, String, String})
precompile(Tuple{typeof(Base.setindex!), NCDatasets.CommonDataModel.CFVariable{Int32, 1, NCDatasets.Variable{Int32, 1, NCDatasets.NCDataset{Nothing, Base.Missing}}, NCDatasets.CommonDataModel.Attributes{NCDatasets.Variable{Int32, 1, NCDatasets.NCDataset{Nothing, Base.Missing}}}, NamedTuple{(:fillvalue, :missing_values, :scale_factor, :add_offset, :calendar, :time_origin, :time_factor, :maskingvalue), Tuple{Nothing, Tuple{}, Nothing, Nothing, Nothing, Nothing, Nothing, Base.Missing}}}, Array{Int32, 1}, Base.OneTo{Int64}})
@static if pkgversion(NCDatasets.DiskArrays) < v"0.4.10"
    precompile(Tuple{typeof(NCDatasets.DiskArrays.writeblock_sizecheck!), NCDatasets.Variable{Int32, 1, NCDatasets.NCDataset{Nothing, Base.Missing}}, Array{Int32, 1}, Base.OneTo{Int64}})
else
    precompile(Tuple{typeof(NCDatasets.DiskArrays.writeblock_checked!), NCDatasets.Variable{Int32, 1, NCDatasets.NCDataset{Nothing, Base.Missing}}, Array{Int32, 1}, Base.OneTo{Int64}})
end
precompile(Tuple{typeof(Base.setindex!), NCDatasets.CommonDataModel.Attributes{NCDatasets.Variable{Float64, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}}, String, String})
precompile(Tuple{typeof(Base.setindex!), NCDatasets.CommonDataModel.CFVariable{Float64, 2, NCDatasets.Variable{Float64, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, NCDatasets.CommonDataModel.Attributes{NCDatasets.Variable{Float64, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}}, NamedTuple{(:fillvalue, :missing_values, :scale_factor, :add_offset, :calendar, :time_origin, :time_factor, :maskingvalue), Tuple{Nothing, Tuple{}, Nothing, Nothing, Nothing, Nothing, Nothing, Base.Missing}}}, Base.ReinterpretArray{Float64, 2, Tuple{Float64, Float64, Float64}, Array{Tuple{Float64, Float64, Float64}, 1}, true}, Base.OneTo{Int64}, Base.OneTo{Int64}})
@static if pkgversion(NCDatasets.DiskArrays) < v"0.4.10"
    precompile(Tuple{typeof(NCDatasets.DiskArrays.writeblock_sizecheck!), NCDatasets.Variable{Float64, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, Base.ReinterpretArray{Float64, 2, Tuple{Float64, Float64, Float64}, Array{Tuple{Float64, Float64, Float64}, 1}, true}, Base.OneTo{Int64}, Vararg{Base.OneTo{Int64}}})
else
    precompile(Tuple{typeof(NCDatasets.DiskArrays.writeblock_checked!), NCDatasets.Variable{Float64, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, Base.ReinterpretArray{Float64, 2, Tuple{Float64, Float64, Float64}, Array{Tuple{Float64, Float64, Float64}, 1}, true}, Base.OneTo{Int64}, Vararg{Base.OneTo{Int64}}})
end
precompile(Tuple{typeof(NCDatasets.DiskArrays.writeblock!), NCDatasets.Variable{Float64, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, Base.ReinterpretArray{Float64, 2, Tuple{Float64, Float64, Float64}, Array{Tuple{Float64, Float64, Float64}, 1}, true}, Base.OneTo{Int64}, Vararg{Base.OneTo{Int64}}})
precompile(Tuple{typeof(NCDatasets._write_data_to_nc), NCDatasets.Variable{Float64, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, Base.ReinterpretArray{Float64, 2, Tuple{Float64, Float64, Float64}, Array{Tuple{Float64, Float64, Float64}, 1}, true}, Base.OneTo{Int64}, Vararg{Base.OneTo{Int64}}})
precompile(Tuple{typeof(Base.to_indices), NCDatasets.Variable{Float64, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, Tuple{Base.OneTo{Int64}, Base.OneTo{Int64}}, Tuple{Base.OneTo{Int64}, Base.OneTo{Int64}}})
precompile(Tuple{typeof(NCDatasets.size_getindex), NCDatasets.Variable{Float64, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, Base.OneTo{Int64}, Vararg{Base.OneTo{Int64}}})
precompile(Tuple{typeof(NCDatasets._size_getindex), NCDatasets.Variable{Float64, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, Tuple{}, Int64, Base.OneTo{Int64}, Base.OneTo{Int64}})
precompile(Tuple{typeof(Base.setindex!), NCDatasets.CommonDataModel.CFVariable{Int32, 2, NCDatasets.Variable{Int32, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, NCDatasets.CommonDataModel.Attributes{NCDatasets.Variable{Int32, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}}, NamedTuple{(:fillvalue, :missing_values, :scale_factor, :add_offset, :calendar, :time_origin, :time_factor, :maskingvalue), Tuple{Nothing, Tuple{}, Nothing, Nothing, Nothing, Nothing, Nothing, Base.Missing}}}, Base.ReinterpretArray{Int32, 2, Tuple{Int32, Int32}, Array{Tuple{Int32, Int32}, 1}, true}, Base.OneTo{Int64}, Base.OneTo{Int64}})
@static if pkgversion(NCDatasets.DiskArrays) < v"0.4.10"
    precompile(Tuple{typeof(NCDatasets.DiskArrays.writeblock_sizecheck!), NCDatasets.Variable{Int32, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, Base.ReinterpretArray{Int32, 2, Tuple{Int32, Int32}, Array{Tuple{Int32, Int32}, 1}, true}, Base.OneTo{Int64}, Vararg{Base.OneTo{Int64}}})
else
    precompile(Tuple{typeof(NCDatasets.DiskArrays.writeblock_checked!), NCDatasets.Variable{Int32, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, Base.ReinterpretArray{Int32, 2, Tuple{Int32, Int32}, Array{Tuple{Int32, Int32}, 1}, true}, Base.OneTo{Int64}, Vararg{Base.OneTo{Int64}}})
end
precompile(Tuple{typeof(NCDatasets.DiskArrays.writeblock!), NCDatasets.Variable{Int32, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, Base.ReinterpretArray{Int32, 2, Tuple{Int32, Int32}, Array{Tuple{Int32, Int32}, 1}, true}, Base.OneTo{Int64}, Vararg{Base.OneTo{Int64}}})
precompile(Tuple{typeof(NCDatasets._write_data_to_nc), NCDatasets.Variable{Int32, 2, NCDatasets.NCDataset{Nothing, Base.Missing}}, Base.ReinterpretArray{Int32, 2, Tuple{Int32, Int32}, Array{Tuple{Int32, Int32}, 1}, true}, Base.OneTo{Int64}, Vararg{Base.OneTo{Int64}}})
for N in 6:10
    precompile(Tuple{typeof(Core.kwcall), NamedTuple{(:force3D, :write_computed), Tuple{Bool, Bool}}, typeof(VoronoiMeshes.save), String, VoronoiMeshes.VoronoiMesh{false, N, Int32, Float64, Zeros.Zero}})
    precompile(Tuple{typeof(Core.kwcall), NamedTuple{(:force3D, :write_computed), Tuple{Bool, Bool}}, typeof(VoronoiMeshes.save), String, VoronoiMeshes.VoronoiMesh{true, N, Int32, Float64, Float64}})

    precompile(Tuple{typeof(VoronoiMeshes.save_to_netcdf!), NCDatasets.NCDataset{Nothing, Base.Missing}, VoronoiMeshes.VoronoiDiagram{false, N, Int32, Float64, Zeros.Zero}, Bool})
    precompile(Tuple{typeof(VoronoiMeshes.save_to_netcdf!), NCDatasets.NCDataset{Nothing, Base.Missing}, VoronoiMeshes.VoronoiDiagram{true, N, Int32, Float64, Float64}, Bool})

    precompile(Tuple{typeof(save_to_netcdf), String, VoronoiMesh{false, N, Int32, Float64, Zeros.Zero}})
    precompile(Tuple{typeof(save_to_netcdf), String, VoronoiMesh{true, N, Int32, Float64, Float64}})
end

