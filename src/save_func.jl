"""
    save(filename::String, mesh::AbstractVoronoiMesh; format=:netcdf5_64bit_data, write_computed=false)

The `NCDatasets` package must be loaded to use this function.
Save the `mesh` in a NetCDF file named `filename`.
By default the CDF5 NetCDF file format is used. Other formats can be specified with the `format` keyword. See other available options with ?[`NCDataset`](@ref).
By default only the base fields are stored in the file. In order to write also the computed fieldin the NetCDF file, use the keywork `write_computed=true`.

---

    save(filename::String, diagram::VoronoiDiagram; format=:netcdf5_64bit_data)

The `NCDatasets` package must be loaded to use this function.
Save the `diagram` in a NetCDF file named `filename`.
By default the CDF5 NetCDF file format is used. Other formats can be specified with the `format` keyword. See other available options with ?[`NCDataset`](@ref).
"""
function save(filename, obj::T; kwds...) where {T<:Union{<:AbstractVoronoiDiagram,<:AbstractVoronoiMesh}}
    name, ext = Base.Filesystem.splitext(filename)
    if ext == ".nc"
        save_to_netcdf(filename, obj; kwds...)
        T<:AbstractVoronoiMesh && write(name*".graph.info", String(take!(graph_partition(obj))))
    else
        error("Unsupported file extension: $filename")
    end
    return nothing
end

for N in 6:9
    precompile(Tuple{typeof(save), String, VoronoiMesh{false, N, Int32, Float64, Zeros.Zero}})
    precompile(Tuple{typeof(save), String, VoronoiMesh{true, N, Int32, Float64, Float64}})
end

