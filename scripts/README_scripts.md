# Scripts

Utility scripts for generating periodic planar Voronoi meshes, computing
per-cell quality metrics, and visualizing both. All scripts are run from this
`scripts/` directory with the local `Project.toml` environment:

```
julia --project=. <script>.jl [args...]
```

All build scripts write into `scripts/output/` (created if missing, and
git-ignored). Filenames are self-descriptive: `mesh_periodic_<kind>_<level>_nc<N>...`.

## Mesh-building scripts

### `build_set_irregular_meshes.jl`

Generates a sequence of same-resolution meshes with increasing irregularity,
starting from a converged centroidal (Level 0) mesh. Only cells within a band
around a fixed line are perturbed at each level; growing Gaussian noise is
added to their generator points, then a few Lloyd iterations are run — enough
to fix degenerate cells without erasing the introduced irregularity. Extra
Lloyd iterations are applied automatically if obtuse Delaunay triangles remain.

```
julia --project=. build_set_irregular_meshes.jl [nc] [num_levels] [base_strength]
julia --project=. build_set_irregular_meshes.jl [nc] 0 [strength]   # single perturbed mesh
```

- `nc` (default 80): target cell count for the Level 0 reference mesh.
- `num_levels` (default 5): number of perturbed levels to build. **Use `0`**
  to build only the Level 0 mesh plus a single perturbation at exactly
  `base_strength` (no multi-level sweep).
- `base_strength` (default 0.15): perturbation amplitude at level 1, as a
  fraction of the average cell spacing dc ≈ 1/√nc. Level `i` uses
  `strength = base_strength * i` (except in the `num_levels=0` case, where
  it's used directly as the single perturbation strength).

Outputs per level: `..._vor.vtu`, `..._tri.vtu`, `....png` (mesh overlay),
`..._<metric>.png` (one per registered metric), plus
`output/irregular_metrics_summary.csv` and `output/irregular_voronoi_meshes.txt`
(manifest, see below).

### `build_set_regular_meshes.jl`

Builds a sequence of regular meshes by successive refinement: Level 0 is a
converged hex mesh; each subsequent level uses the previous level's cell
centers + edge midpoints as generators and re-converges via Lloyd relaxation,
roughly quadrupling the cell count each time.

```
julia --project=. build_set_regular_meshes.jl [nc_ref] [num_levels]
```

- `nc_ref` (default 16): target cell count for the Level 0 mesh.
- `num_levels` (default 4): number of refinement levels.

Same outputs as above, under `output/regular_metrics_summary.csv` /
`output/regular_voronoi_meshes.txt`.

### `build_set_refined_meshes_vtu.jl`

Builds a sequence of independently-generated centroidal meshes with cell
counts growing as `base_cells * 2^i` (geometric refinement, not derived from
one another).

```
julia --project=. build_set_refined_meshes_vtu.jl [base_cells] [num_scales] [ini_scale]
```

- `base_cells` (default 16), `num_scales` (default 11), `ini_scale` (default 0):
  cell counts are `base_cells*2^ini_scale, ..., base_cells*2^(ini_scale+num_scales-1)`
  (defaults give 16, 32, ..., 16384).

Same per-mesh outputs as above, under `output/refined_metrics_summary.csv` /
`output/refined_voronoi_meshes.txt`.

All three build scripts above share the same structure: a `main(...)`
function driven by positional CLI args (with the same defaults as running
with no args), building each level/scale through the shared helpers below,
then a single call to `MeshTools.finalize_mesh_set` to print/save the summary
table, CSV, and manifest.

## Metrics and plotting

### `mesh_tools.jl`

Not run directly — the shared build/report pipeline, `include`d by the three
build scripts above and by `plot_mesh_properties.jl`. Defines:

- Per-cell metric functions (each mesh → one `Vector` per cell), registered in
  `METRICS`: `cell_area_normalized`, `cell_distortion`, `cell_distortion_rms`,
  `cell_diameter_normalized`, `cell_alignment`. Add a new metric by writing
  the function and adding it to the `METRICS` tuple — it then automatically
  shows up everywhere `METRICS` is iterated (printing, plotting, CSV export).
- Shared plotting helpers (need a Makie backend loaded first, e.g.
  `using GLMakie`): `save_mesh_png` (plain Voronoi + triangulation overlay),
  `save_property_png` / `save_all_metric_pngs` (per-cell-metric colored plots).
- Shared printing/export helpers: `print_metrics_summary`, `metrics_summary_row`,
  `print_summary_table`, `save_summary_csv`, `save_manifest`, `rebuild_manifest`
  (scans an output directory for every mesh file matching a given kind's
  naming pattern and writes the manifest from that — not just the meshes the
  current run produced — so a manifest stays complete even when it should
  include meshes left over from an earlier run).
- Shared build-script helpers, factored out of what used to be duplicated
  across the three build scripts: `build_hex_reference` (the common Level-0
  regular hex mesh, rebuilt against the exact periodic domain),
  `save_mesh_level` (saves one mesh's VTU + overlay PNG + per-metric PNGs +
  metrics summary, returning its summary row), and `finalize_mesh_set` (the
  common end-of-build step: summary table + CSV + manifest rebuild).

### `plot_mesh_properties.jl`

Reads previously-saved meshes (VTU) and, for each, computes every metric in
`MeshTools.METRICS`, plots it as a colored Voronoi diagram, and prints/saves
a combined summary table.

```
julia --project=. plot_mesh_properties.jl                     # process every manifest in output/
julia --project=. plot_mesh_properties.jl output/some_manifest.txt
julia --project=. plot_mesh_properties.jl output/mesh_vor.vtu  # single mesh
```

With no argument, it looks for `output/irregular_voronoi_meshes.txt`,
`output/regular_voronoi_meshes.txt`, and `output/refined_voronoi_meshes.txt` —
the manifests written automatically by the three build scripts above (each
lists one `*_vor.vtu` filename per line, relative to the manifest's own
directory). This is the easiest way to (re)generate metric plots for
everything currently in `output/` without hunting down individual files.

Output: `output/<mesh>_<metric>.png` for every mesh/metric pair, plus a
printed table and `output/mesh_properties_summary.csv`.

## Other scripts

- `save_regular_mesh_vtu.jl` — minimal worked example of building a small
  mesh, saving/reading Voronoi + triangulation VTU files, and plotting.
  Not part of the metrics workflow above.
- `create_distorted_meshes.jl` — older standalone script that builds a
  sequence of variable-resolution meshes (density-function-driven Lloyd
  relaxation) and saves them to NetCDF. Usage:
  `julia --project=. create_distorted_meshes.jl <x_period> <y_period> <dc>`.
  Predates `mesh_tools.jl` and does not use it.

## Typical workflow

```
julia --project=. build_set_regular_meshes.jl
julia --project=. build_set_irregular_meshes.jl
julia --project=. build_set_refined_meshes_vtu.jl
julia --project=. plot_mesh_properties.jl        # metric plots + summary CSV for everything above
```
