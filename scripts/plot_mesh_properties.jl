# plot_mesh_properties.jl
#
# Reads one or more previously-saved Voronoi meshes (VTU format) and, for each
# one, computes per-cell area, distortion and diameter, then plots each
# property as a colored Voronoi diagram saved to a PNG.
#
# Usage:
#   julia --project=. plot_mesh_properties.jl [manifest.txt | mesh_vor.vtu]
#
#   With no argument, processes every "output/*_voronoi_meshes.txt" manifest
#   produced by the build_set_*.jl scripts.
#
#   - manifest.txt : a text file listing one "*_vor.vtu" path per line
#                    (blank lines and lines starting with '#' are ignored;
#                    relative paths are resolved against the manifest's
#                    own directory). This is the format written by the
#                    build_set_*.jl scripts into scripts/output/*.txt.
#   - mesh_vor.vtu : a single mesh, given either as "<base>_vor.vtu" (the
#                    matching "<base>_tri.vtu" is derived automatically) or
#                    as a bare "<base>.vtu" base name.
#
# Output: for each mesh "<base>", writes "<base>_area.png",
# "<base>_distortion.png" and "<base>_diameter.png" next to the input file.
#
# Note: PNG export uses GLMakie. On headless servers swap to CairoMakie
# (add it to the environment and replace `using GLMakie` below).

using VoronoiMeshes
using TensorsLite
using TensorsLiteGeometry
using ReadVTK
using GLMakie
using Statistics

# Per-cell area, already computed by the package.
cell_area(mesh) = mesh.cells.area

# Per-cell irregularity: (max_edge - min_edge) / mean_edge over each cell's edges.
function cell_distortion(mesh)
    edge_lengths    = mesh.edges.length
    cells_edges     = mesh.cells.edges
    nEdges_per_cell = mesh.cells.nEdges
    nc = length(nEdges_per_cell)

    d = Vector{Float64}(undef, nc)
    for c in 1:nc
        ne   = Int(nEdges_per_cell[c])
        lmin = Inf; lmax = -Inf; lsum = 0.0
        for j in 1:ne
            l    = edge_lengths[cells_edges[c][j]]
            lmin = min(lmin, l)
            lmax = max(lmax, l)
            lsum += l
        end
        d[c] = (lmax - lmin) / (lsum / ne)
    end
    return d
end

# Per-cell diameter: max distance between any two of the cell's vertices,
# using their periodic images closest to the cell center.
function cell_diameter(mesh)
    vert_pos     = mesh.vertices.position
    cell_pos     = mesh.cells.position
    verticesOnCell = mesh.cells.vertices
    nEdges_per_cell = mesh.cells.nEdges
    xp, yp = mesh.x_period, mesh.y_period
    nc = length(cell_pos)

    diam = Vector{Float64}(undef, nc)
    for c in 1:nc
        ne = Int(nEdges_per_cell[c])
        cp = cell_pos[c]
        vs = verticesOnCell[c]
        local_vertices = ntuple(j -> closest(cp, vert_pos[vs[j]], xp, yp), ne)
        dmax = 0.0
        for i in 1:ne, j in (i + 1):ne
            dv = local_vertices[i] - local_vertices[j]
            dmax = max(dmax, sqrt(dv.x^2 + dv.y^2))
        end
        diam[c] = dmax
    end
    return diam
end

function load_mesh(path)
    if endswith(path, "_vor.vtu")
        tri_path = replace(path, "_vor.vtu" => "_tri.vtu")
        return VoronoiMeshes.read_from_vtu(path, tri_path)
    else
        return VoronoiMeshes.read_from_vtu(path)
    end
end

function base_label(path)
    name = basename(path)
    name = replace(name, "_vor.vtu" => "", ".vtu" => "")
    return joinpath(dirname(path), name)
end

function save_property_png(filename, mesh, values, title)
    polygons = VoronoiMeshes.create_cell_polygons(mesh)
    fig = Figure(size = (750, 700))
    ax  = Axis(fig[1, 1], title = title, aspect = DataAspect())
    plt = poly!(ax, polygons, color = values, colormap = :viridis, strokewidth = 0.5, strokecolor = (:black, 0.3))
    Colorbar(fig[1, 2], plt)
    hidedecorations!(ax)
    GLMakie.save(filename, fig)
    return nothing
end

function process_mesh(path)
    println("Reading: $path")
    mesh = load_mesh(path)
    nc = length(mesh.cells.position)
    println("  nc = $nc")

    properties = (
        ("area", cell_area(mesh)),
        ("distortion", cell_distortion(mesh)),
        ("diameter", cell_diameter(mesh)),
    )

    label = base_label(path)
    for (name, values) in properties
        println("  $name: mean = $(round(mean(values), digits=5)), min = $(round(minimum(values), digits=5)), max = $(round(maximum(values), digits=5))")
        filename = "$(label)_$(name).png"
        save_property_png(filename, mesh, values, "$(basename(label)) — $name")
        println("  Saved: $filename")
    end
    return nothing
end

function resolve_mesh_paths(input_path)
    if endswith(input_path, ".txt")
        dir = dirname(input_path)
        lines = readlines(input_path)
        paths = String[]
        for line in lines
            l = strip(line)
            (isempty(l) || startswith(l, "#")) && continue
            push!(paths, isabspath(l) ? l : joinpath(dir, l))
        end
        return paths
    else
        return [input_path]
    end
end

# Default: no argument given -> process every manifest already in output/.
default_manifests() = filter(isfile, joinpath.("output", (
    "regular_voronoi_meshes.txt", "irregular_voronoi_meshes.txt", "refined_voronoi_meshes.txt",
)))

input_paths = if length(ARGS) >= 1
    [ARGS[1]]
else
    manifests = default_manifests()
    isempty(manifests) && error("Usage: julia --project=. plot_mesh_properties.jl <manifest.txt | mesh_vor.vtu>")
    println("No argument given, defaulting to manifests: $(join(manifests, ", "))")
    manifests
end

for input_path in input_paths, path in resolve_mesh_paths(input_path)
    process_mesh(path)
end

println("\nDone!")
