# Scripts

Utility scripts for generating periodic planar Voronoi meshes and their
per-cell quality metrics. Run from this `scripts/` directory:

```
julia --project=. <script>.jl [args...]
```

Every `build_set_*.jl` script writes to its own `output/<kind>_<params>/`
subfolder and accepts `--help`/`-h` for full usage.

## Mesh-building scripts

### `build_set_global_distortion_meshes.jl`

Same-resolution meshes with increasing distortion: every cell's generator
point is randomly perturbed, then lightly relaxed to remove bad triangles.

```
julia --project=. build_set_global_distortion_meshes.jl [nc] [num_levels] [base_strength] [fixup_iters]
```

- `nc` (default 80): cell count, fixed across all levels.
- `num_levels` (default 6): number of perturbed levels.
- `base_strength` (default 0.05): perturbation amplitude at level 1, as a
  fraction of the mean cell spacing; level `i` uses `base_strength * i`.
- `fixup_iters` (default 4): Lloyd iterations applied after perturbing, to
  clean up obtuse triangles. Raise this if obtuse triangles remain a problem
  at the strengths you need.

### `build_set_regular_meshes.jl`

Regular meshes by successive refinement, quadrupling cell count each level.

```
julia --project=. build_set_regular_meshes.jl [nc_ref] [num_levels]
```

- `nc_ref` (default 64): cell count for Level 0.
- `num_levels` (default 4): number of levels (Level 0 + refinements) — same
  convention as `circular_refined`'s `num_levels` and `refined_meshes_vtu`'s
  `num_scales`, so the same number gives the same finest cell count on all three.

### `build_set_refined_meshes_vtu.jl`

Independently-generated meshes with cell counts growing as `base_cells * 4^i`
(not derived from one another) — same cell-count ladder as `regular`/`circular_refined`.

```
julia --project=. build_set_refined_meshes_vtu.jl [base_cells] [num_scales] [ini_scale]
```

- `base_cells` (default 64), `num_scales` (default 4), `ini_scale` (default 0).

### `build_set_circular_refined_meshes.jl`

Locally-refined meshes: cell density is 16x higher within an off-center
circular region than in the background, tapering smoothly. Refines by
quadrupling cell count each level, like `build_set_regular_meshes.jl`.

```
julia --project=. build_set_circular_refined_meshes.jl [nc_ref] [num_levels]
```

- `nc_ref` (default 64): cell count for Level 0.
- `num_levels` (default 4): number of levels (Level 0 + refinements).

All four scripts write per level: `..._vor.vtu`, `..._tri.vtu`, `....png`
(mesh overlay), `..._<metric>.png` (one per metric), a `<kind>_metrics_summary.csv`,
and a `<kind>_voronoi_meshes.txt` manifest.

## Metrics and plotting

### `mesh_tools.jl`

Shared library (metric functions, plotting, build/report helpers) used by
all scripts above. Not run directly.

### `plot_mesh_properties.jl`

Recomputes metrics and per-metric PNGs for already-saved meshes.

```
julia --project=. plot_mesh_properties.jl                     # every manifest under output/
julia --project=. plot_mesh_properties.jl output/some_manifest.txt
julia --project=. plot_mesh_properties.jl output/mesh_vor.vtu  # single mesh
```

### `plot_metrics_summary.jl`

Plots a convergence figure from one `<kind>_metrics_summary.csv` — cell
count on the x-axis if it varies row to row, else perturbation strength `d`.

```
julia --project=. plot_metrics_summary.jl output/<run>/<kind>_metrics_summary.csv
```

Output: `<csv_basename>_convergence.pdf`/`.eps` next to the input CSV.

## Other scripts

- `save_regular_mesh_vtu.jl` — minimal worked example (build, save/read VTU, plot).
- `create_distorted_meshes.jl` — older standalone script, predates `mesh_tools.jl`.
  Usage: `julia --project=. create_distorted_meshes.jl <x_period> <y_period> <dc>`.

## Typical workflow

```
julia --project=. build_set_regular_meshes.jl
julia --project=. build_set_global_distortion_meshes.jl
julia --project=. build_set_refined_meshes_vtu.jl
julia --project=. build_set_circular_refined_meshes.jl
julia --project=. plot_mesh_properties.jl        # metric plots + summary CSV for everything above
```

## Production runs (cluster)

Convergence cases share the same 64→65536 cell-count ladder
(64/256/1024/4096/16384/65536), distortion case at 10k cells, all run in
parallel:

```bash
julia -O3 --threads=2 --project=. build_set_regular_meshes.jl 64 6               > regular.log 2>&1 &
julia -O3 --threads=2 --project=. build_set_refined_meshes_vtu.jl 64 6           > refined.log 2>&1 &
julia -O3 --threads=2 --project=. build_set_circular_refined_meshes.jl 64 6      > circular_refined.log 2>&1 &
julia -O3 --threads=2 --project=. build_set_global_distortion_meshes.jl 10000 6 0.05 10  > global_distortion.log 2>&1 &

wait
```

Finest cell counts: `regular`/`refined`/`circular_refined` all =65536,
`global_distortion` =10000 (fixed resolution, distortion sweep).
`circular_refined` is dramatically slower than the others (density-weighted
Lloyd) — expect it to still be running after the rest finish (observed
~2.5h for this exact command on a desktop machine).

`global_distortion`'s `fixup_iters=10` (vs. the script's own default of 4)
was validated overnight: it keeps obtuse triangles near zero across the
whole sweep, at the cost of roughly halving the achieved distortion at each
nominal `d` (more Lloyd relaxation smooths the perturbation back down). If
you need comparable distortion magnitude *and* clean triangles, raise
`base_strength` alongside `fixup_iters` rather than using `fi=10` with the
default strengths.

