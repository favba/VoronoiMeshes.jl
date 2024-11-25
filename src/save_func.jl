is_netcdf_ext(s::AbstractString) = last(s,3) == ".nc"

function save(filename, obj::T; kwds...) where {T<:Union{<:AbstractVoronoiDiagram,<:VoronoiMesh}}
    if is_netcdf_ext(filename)
        save_to_netcdf(filename, obj; kwds...)
    else
        error("Unsupported file extension: $filename")
    end
    return nothing
end

