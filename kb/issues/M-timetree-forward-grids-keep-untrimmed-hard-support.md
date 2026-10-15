# Forward-pass distributions keep their full hard-side support

Forward-pass grids hold a median 17–124 times more points than their mass core, all on the hard left side, where the density is below $10^{-50}$ of the peak. Forward convolution, division and multiplication take about 95% of the distribution-operation time and about 55% of the wall time of a 400-tip run. Measurements and mechanism: [kb/proposals/timetree-forward-pass-grid-support.md](../proposals/timetree-forward-pass-grid-support.md).

## Locations

- `fn lower_edge()`, `fn upper_edge()` in [packages/treetime-distribution/src/distribution_ops/mass_domain.rs](../../packages/treetime-distribution/src/distribution_ops/mass_domain.rs): a `Hard` or `HardApproach` side returns the grid edge, regardless of mass
- `fn convolve_across_edge()` in [packages/treetime-distribution/src/distribution_ops/edge_convolution.rs](../../packages/treetime-distribution/src/distribution_ops/edge_convolution.rs): the window is clipped to the result's own edge on the hard side
- `fn refine_distribution_from_parent()` in [packages/treetime/src/timetree/inference/forward_pass.rs](../../packages/treetime/src/timetree/inference/forward_pass.rs): the posterior (forward message × subtree distribution) spans the full support intersection and is not re-windowed

## Fix direction

Trim both sides of the forward convolution output and of the forward posterior at $\epsilon$ of the mass, hard sides included. Projected median forward FFT length: 34 676 → 719 (mpox 1000), 85 626 → 1 971 (400-tip simulation).

> [!IMPORTANT]
> **Decision required.** What a trimmed hard side becomes. Options: (a) a `Linear` tail fitted at the new edge, which bridges gaps under the empty invariant but drops the hard bound; (b) a `Linear` tail truncated at the original bound, which keeps the exact support but needs a new `BoundaryBehavior` variant handled by multiplication, division, mass and tail composition; (c) `Hard` at the trimmed edge, which moves a bound inward. Evidence: shifting a hard edge inward by about one narrow-operand mean produced an empty posterior on `data/sc2/4500`. Options (a) and (b) change the tail assignments in [kb/decisions/distribution-tails-and-arithmetic.md](../decisions/distribution-tails-and-arithmetic.md) (forward message: `Hard` left).

## Validation

- `cargo test` including the golden-master runner tests
- Validation suite (`validation/`, `--suite timetree`): date RMSE per (L, mu) unchanged against base; no failures; no negative branches
- `data/mpox/clade-ii/1000`, `data/sc2/4500`, `data/ebola/362`, `data/rsv/a/100`: no node without a date, no child dated before its parent, node dates against base
- Per-operation grid sizes and time before and after, with the measurement set-up of the proposal
