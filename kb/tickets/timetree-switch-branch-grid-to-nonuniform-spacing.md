# Switch branch distribution grid to non-uniform spacing

The branch distribution grid in `create_simple_grid()` uses `Array1::linspace` (uniform spacing) over the full range `[min_bl, max(peak_max_bl, MAX_BRANCH_TIME * clock_rate)]`. For short branches where `MAX_BRANCH_TIME * clock_rate` dominates, the peak region receives fewer grid points than needed for high-accuracy time estimates.

## Quantitative impact

For a 1-mutation flu branch (center=3.3e-4, rate=0.003, L=3000):

- Current grid (1000 points, wide range): dx=6e-4, ~0.56 points per sigma
- Previous grid (200 points, narrow range): dx=1.65e-5, ~20 points per sigma

The argmax error is bounded by dx/2 = 3e-4 subs/site = 0.1 years per branch. Errors accumulate through the tree: golden master tests comparing v1 against v0 require 1.0-year tolerance, with the worst-case node (root-adjacent) differing by 0.92 years. Validation RMSE across all nodes is 0.091 years.

Peak shape is under-resolved for confidence interval estimation.

## Affected tests

Golden master tests `test_gm_runner_marginal_dense` and `test_gm_runner_marginal_sparse` use `epsilon = 1e0` (1.0 years). With a non-uniform grid restoring peak resolution, this tolerance should tighten.

## v0 reference

v0 uses a 5-segment non-uniform grid (`branch_len_interpolator.py:50-62`):

- Log-spaced near zero (5 points)
- Linear near zero (8 points)
- Quadratic from 0 to peak (n/3 points, dense near peak)
- Quadratic from peak to 3\*sigma (n/3 points)
- Quadratic from 3\*sigma to MAX_BRANCH_LENGTH (n/3 points, sparse tail)

This concentrates ~40 of 125 points near the peak while covering up to MAX_BRANCH_LENGTH = 4.0 subs/site.

## Fix

Implement a non-uniform grid matching v0's multi-segment design, or use a two-tier approach: fine uniform grid over the peak region concatenated with a coarser grid for the tail. Requires `DistributionFunction` to support non-uniform grids or a resample step.

## Related issues

- Source: [kb/issues/M-timetree-branch-grid-uniform-resolution.md](../issues/M-timetree-branch-grid-uniform-resolution.md) -- delete after full resolution
