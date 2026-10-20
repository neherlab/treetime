# Timetree grid spacing never coarsens below the finest operand

`fn resample_to_mass_window()` [packages/treetime-distribution/src/distribution_ops/mass_domain.rs#L97-L125](../../packages/treetime-distribution/src/distribution_ops/mass_domain.rs#L97-L125) uses the smaller of the window spacing (window / (`GRID_POINTS` - 1)) and the spacing of its input. In `fn convolve_across_edge()` [packages/treetime-distribution/src/distribution_ops/edge_convolution.rs#L10-L43](../../packages/treetime-distribution/src/distribution_ops/edge_convolution.rs#L10-L43) the input is the FFT result on the finer operand's spacing, and multiplication results are never coarser than the finest operand (grid intersection contract, [kb/decisions/distribution-tails-and-arithmetic.md](../decisions/distribution-tails-and-arithmetic.md#grid-intersection-contract)). The finest spacing of any factor therefore propagates to every later distribution.

On base `94f921f8e7`, mass-core points over `GRID_POINTS` = 300, median / 90th percentile of forward outputs: mpox 1000 0.6-2.0 / 1.3-3.0, 400-tip simulation 2.3-5.0 / 3.4-8.4. After trimming to mass cores, operands on the coarser instead of the finer spacing shorten the median forward FFT from 719 to 452 points (mpox) and from 1 971 to 417 (simulation). Measurements: [kb/proposals/timetree-forward-pass-grid-support.md](../proposals/timetree-forward-pass-grid-support.md).

## Fix direction

Cap the convolution output spacing at the coarser operand's spacing instead of the FFT result's: `dx = min(window / (GRID_POINTS - 1), max(dx_a, dx_b))`.

The experiment commit `6d06d211a` on branch `claude-grid-trim` prototypes the cap behind the `TT_COARSE` environment variable (`max_dx` argument of its `fn resample_trimmed()`), together with the forward-posterior trimming of [M-timetree-forward-grids-keep-untrimmed-hard-support.md](M-timetree-forward-grids-keep-untrimmed-hard-support.md). With both changes, every node of `data/ebola/362`, `data/mpox/clade-ii/1000` and `data/sc2/4500` keeps a date, the largest node-date change against base is 0.027 years (mpox), and the 2880-run validation suite has no failures. Results: [kb/proposals/timetree-forward-pass-grid-support.md](../proposals/timetree-forward-pass-grid-support.md#experiment-results-2026-10-19).

> [!IMPORTANT]
> **Investigation required.** An earlier implementation of the cap, on the removed branch `moment-matched-convolution` (`max_dx` argument of `fn resample_to_mass_window()`), moved 60 node dates of `data/mpox/clade-ii/1000` by more than a year against base. That code no longer exists, so the result cannot be reproduced, and the prototype above does not show it. Before the cap lands, compare the node dates of `data/mpox/clade-ii/1000` with and without the cap, once with and once without the posterior trimming, to find out whether the cap alone causes such shifts.

Land the cap in the same change as the posterior trimming of [M-timetree-forward-grids-keep-untrimmed-hard-support.md](M-timetree-forward-grids-keep-untrimmed-hard-support.md). Without the cap, a trimmed posterior gridded at 300 points over its mass window passes its fine spacing into the next convolution with a wide branch distribution, and the forward grids grow (sc2: median forward convolution 9 989 to 13 349 points, maximum 240 050; mpox: a grid of 1.5 million points).

## Validation

As in [M-timetree-forward-grids-keep-untrimmed-hard-support.md](M-timetree-forward-grids-keep-untrimmed-hard-support.md); on the validation suite, compare date RMSE and confidence-interval coverage per (L, mu) against base.
