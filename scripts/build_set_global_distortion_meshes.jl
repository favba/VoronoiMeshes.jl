using VoronoiMeshes
using DelaunayTriangulation
using TensorsLite
using WriteVTK
using GLMakie

include("mesh_tools.jl")
using .MeshTools

const X_PERIOD = 1.0
const Y_PERIOD = 1.0

# Clamped in absolute terms (not scaled by strength) so a displacement can
# never get large enough to invert/overlap a cell.
const MAX_DISPLACEMENT_FRAC = 0.3

const DEFAULT_FIXUP_ITERS = 4

const MESH_PATTERN = r"^mesh_periodic_global_distortion_nc\d+_d([\d.]+)_vor\.vtu$"

function print_obtuse_triangles(mesh, label)
    n_obtuse, n_tri = MeshTools.obtuse_triangle_count(mesh)
    println("    $label: obtuse triangles = $n_obtuse / $n_tri")
    return n_obtuse
end

function perturb_positions(mesh, strength)
    pos = mesh.cells.position
    nc = length(pos)
    dc = sqrt(X_PERIOD * Y_PERIOD / nc)
    max_disp = MAX_DISPLACEMENT_FRAC * dc

    x_new = copy(pos.x)
    y_new = copy(pos.y)
    for c in 1:nc
        p = pos[c]
        dx = clamp(strength * dc * randn(), -max_disp, max_disp)
        dy = clamp(strength * dc * randn(), -max_disp, max_disp)
        x_new[c] = mod(p.x + dx, X_PERIOD)
        y_new[c] = mod(p.y + dy, Y_PERIOD)
    end
    return VecArray(x = x_new, y = y_new)
end

function build_reference_mesh(nc, outdir)
    mesh, dc = MeshTools.build_hex_reference(nc, X_PERIOD, Y_PERIOD)
    actual_nc = length(mesh.cells.position)
    println("Level 0: creating regular hex mesh (reference nc≈$nc, actual=$actual_nc, dc=$(round(dc, digits=4)))...")

    label = "mesh_periodic_global_distortion_nc$(actual_nc)_d0.0"
    row = MeshTools.save_mesh_level(outdir, mesh, label)
    print_obtuse_triangles(mesh, "Level 0")

    return mesh, actual_nc, row
end

# Always perturbs `ref_mesh` (Level 0), not the previous level's result, so
# each level is an independent draw at its own strength rather than a
# compounding chain. fixup_iters is applied unconditionally (not only when
# obtuse triangles appear) so every level gets the same treatment —
# conditional fixup made the distortion sweep non-monotonic, since a level
# that happened to need no cleanup kept more of its raw distortion than
# neighboring levels that did.
function perturb_level(ref_mesh, strength, level, actual_nc, outdir, fixup_iters)
    println("Level $level: perturbing (strength=$(round(strength, digits=3)))...")

    generators = perturb_positions(ref_mesh, strength)
    raw_mesh = VoronoiMesh(generators, X_PERIOD, Y_PERIOD; max_iter = 0)
    print_obtuse_triangles(raw_mesh, "before fixup")
    mesh = VoronoiMesh(generators, X_PERIOD, Y_PERIOD; max_iter = fixup_iters)
    print_obtuse_triangles(mesh, "after fixup")

    label = "mesh_periodic_global_distortion_nc$(actual_nc)_d$(round(strength, digits=3))"
    return MeshTools.save_mesh_level(outdir, mesh, label)
end

function main(nc, num_levels, base_strength, fixup_iters)
    outdir = MeshTools.run_outdir("global_distortion_nc$(nc)_L$(num_levels)_bs$(round(base_strength, digits=3))_fi$(fixup_iters)")
    MeshTools.save_run_info(outdir, [
        "Script: build_set_global_distortion_meshes.jl",
        "",
        "CLI parameters:",
        "  nc = $nc",
        "  num_levels = $num_levels",
        "  base_strength = $base_strength",
        "  fixup_iters = $fixup_iters",
        "",
        "Domain: periodic $(X_PERIOD) x $(Y_PERIOD)",
        "Density: uniform (no density function) for the Level 0 reference mesh.",
        "",
        "Distortion mechanism: every cell's generator point is independently displaced",
        "(Gaussian, sigma = strength*dc, clamped to +/-$(MAX_DISPLACEMENT_FRAC)*dc) from Level 0's",
        "original positions each level. strength = base_strength * level.",
        "$(fixup_iters) Lloyd iterations are always applied afterward.",
    ])

    mesh, actual_nc, row0 = build_reference_mesh(nc, outdir)
    rows = [row0]

    if num_levels == 0
        push!(rows, perturb_level(mesh, base_strength, 1, actual_nc, outdir, fixup_iters))
    else
        for i in 1:num_levels
            strength = base_strength * i
            push!(rows, perturb_level(mesh, strength, i, actual_nc, outdir, fixup_iters))
        end
    end

    MeshTools.finalize_mesh_set(outdir, "global_distortion", rows, MESH_PATTERN, MeshTools.numeric_sort_key(MESH_PATTERN))
    return nothing
end

const USAGE = """
Usage: julia --project=. build_set_global_distortion_meshes.jl [nc] [num_levels] [base_strength] [fixup_iters]
       julia --project=. build_set_global_distortion_meshes.jl [nc] 0 [strength]   # single perturbed mesh

Generates a series of same-resolution meshes with increasing distortion, by
randomly perturbing every cell's generator point and lightly relaxing away
any resulting obtuse Delaunay triangles.

Arguments (all optional, positional):
  nc             Target cell count for the Level 0 reference mesh (default 80).
  num_levels     Number of perturbed levels to build (default 6). Use 0 to
                 build only Level 0 plus a single perturbation at exactly
                 base_strength (no multi-level sweep).
  base_strength  Perturbation amplitude at level 1, as a fraction of the
                 average cell spacing dc ≈ 1/√nc (default 0.05). Level i uses
                 strength = base_strength * i (default sweep: d=0.05,0.1,...,0.3).
  fixup_iters    Lloyd iterations always applied after perturbing, to clean
                 up obtuse triangles (default 4).
"""
MeshTools.handle_help(ARGS, USAGE)

nc = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 80
num_levels = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 6
base_strength = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 0.05
fixup_iters = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : DEFAULT_FIXUP_ITERS
base_strength < 0 && error("base_strength must be >= 0 (got $base_strength)")

main(nc, num_levels, base_strength, fixup_iters)
