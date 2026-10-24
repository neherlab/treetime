# Forward messages that resolve their posterior can build grids far larger than the untrimmed run

The forward pass computes the message from the parent on a capped spacing and checks the posterior it forms. When that posterior holds fewer than `GRID_POINTS - 1` cells across its mass window, `fn ForwardConvolution::resolve_posterior()` resamples the message once at a finer spacing, and then at the fine spacing of the convolution when the posterior is still under-resolved ([packages/treetime-distribution/src/distribution_ops/edge_convolution.rs#L51-L114](../../packages/treetime-distribution/src/distribution_ops/edge_convolution.rs#L51-L114)). Each stage resamples the whole convolution window at the chosen spacing through `fn resample_to_mass_window()` ([packages/treetime-distribution/src/distribution_ops/mass_domain.rs#L332-L368](../../packages/treetime-distribution/src/distribution_ops/mass_domain.rs#L332-L368)), so the point count is the window width divided by that spacing, bounded only by the grid point limit `MaxGridPoints` (default 1 000 000, [packages/treetime-grid/src/max_grid_points.rs#L20](../../packages/treetime-grid/src/max_grid_points.rs#L20)).

## Evidence

- `data/mpox/clade-ii/1000` at `--clock-rate=5e-4` (2026-10-24): at "Resolving spacing" the forward message of leaf `OP696840` resamples its whole 11.6-year window to 43 513 points, while its refined posterior has 301 points
- An earlier example, `data/flu/h3n2/200` with `--clock-rate=3.058e-3 --keep-root` (67 449 points), is not evidence: that input tree is rooted on a 2007 tip, so `--keep-root` contradicts the 1968 to 1972 samples
- The A6 measurements of [kb/decisions/timetree-forward-posterior-trimming.md](../decisions/timetree-forward-posterior-trimming.md) show no run time or memory cost of this size on the default runs of `data/sc2/4500`, `data/ebola/362`, `data/flu/h3n2/200` and the 36 validation runs
- No run has reached the grid point limit through this path

> [!IMPORTANT]
> **Investigation required.** Which stage builds the 43 513-point grid of the mpox case (the resample at the posterior's resolution or the fine fallback), and why the window it resamples is so much wider than the windows of the untrimmed run, for example a window that stays at the full convolution support while the posterior is narrow. Record the stage counts and the largest forward grid per dataset across the smoke datasets and the validation runs, and check whether any input brings the forward pass near `MaxGridPoints`. A finding of that kind needs a bound on the resampled window, for example the posterior's mass window with a margin, instead of the full convolution window.

## Related

- [M-timetree-narrow-branch-distribution-forces-millions-of-grid-points.md](M-timetree-narrow-branch-distribution-forces-millions-of-grid-points.md): another path to grids near the limit
