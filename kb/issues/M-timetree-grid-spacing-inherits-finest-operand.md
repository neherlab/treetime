# Timetree grid spacing never coarsens below the finest operand

`fn resample_to_mass_window()` [packages/treetime-distribution/src/distribution_ops/mass_domain.rs#L107-L121](../../packages/treetime-distribution/src/distribution_ops/mass_domain.rs#L107-L121) uses the smaller of the window spacing (window / (`GRID_POINTS` - 1)) and the spacing of its input. In `fn convolve_across_edge()` [packages/treetime-distribution/src/distribution_ops/edge_convolution.rs#L18-L42](../../packages/treetime-distribution/src/distribution_ops/edge_convolution.rs#L18-L42) the input is the FFT result on the finer operand's spacing, and multiplication results are never coarser than the finest operand (grid intersection contract, [kb/decisions/distribution-tails-and-arithmetic.md](../decisions/distribution-tails-and-arithmetic.md#grid-intersection-contract)). The finest spacing of any factor therefore propagates to every later distribution.

On base `94f921f8e7`, mass-core points over `GRID_POINTS` = 300, median / 90th percentile of forward outputs: mpox 1000 0.6-2.0 / 1.3-3.0, 400-tip simulation 2.3-5.0 / 3.4-8.4. After trimming to mass cores, operands on the coarser instead of the finer spacing shorten the median forward FFT from 719 to 452 points (mpox) and from 1 971 to 417 (simulation). Measurements: [kb/proposals/timetree-forward-pass-grid-support.md](../proposals/timetree-forward-pass-grid-support.md).

## Fix direction

Cap the convolution output spacing at the coarser operand's spacing instead of the FFT result's: `dx = min(window / (GRID_POINTS - 1), max(dx_a, dx_b))`. Implemented and tested on the unmerged branch `moment-matched-convolution` (`max_dx` argument of `fn resample_to_mass_window()`); on `data/mpox/clade-ii/1000` it moved 60 node dates by more than a year against base, with no ground truth to decide which is right.

Land after [M-timetree-forward-grids-keep-untrimmed-hard-support.md](M-timetree-forward-grids-keep-untrimmed-hard-support.md) and measure separately.

## Validation

As in [M-timetree-forward-grids-keep-untrimmed-hard-support.md](M-timetree-forward-grids-keep-untrimmed-hard-support.md); on the validation suite, compare date RMSE and confidence-interval coverage per (L, mu) against the trimmed-support baseline.
