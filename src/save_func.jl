is_netcdf_ext(s::AbstractString) = last(s,3) == ".nc"

"""
    save(filename::String, mesh::VoronoiMesh; format=:netcdf5_64bit_data, write_computed=false)

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
function save(filename, obj::T; kwds...) where {T<:Union{<:AbstractVoronoiDiagram,<:VoronoiMesh}}
    if is_netcdf_ext(filename)
        save_to_netcdf(filename, obj; kwds...)
    else
        error("Unsupported file extension: $filename")
    end
    return nothing
end

