# Switch branch distribution grid to non-uniform spacing

The branch distribution grid in `create_simple_grid()` uses `Array1::linspace` (uniform spacing) over `[one_mutation * 0.01, min(center * 5, MAX_BRANCH_LENGTH)]` with 300 points. For short branches where `one_mutation * 10` governs the extent, the uniform grid may under-resolve the peak compared to v0's non-uniform grid.

## Acceptance criterion

Un-ignore the golden master tests `test_gm_runner_marginal_dense`, `test_gm_runner_marginal_sparse`, and `test_gm_runner_poisson` and verify they pass at `epsilon = 1e-6`. If they do, the uniform grid is sufficient and this ticket can close without a non-uniform implementation.

## v0 reference

v0 uses a 5-segment non-uniform grid (`branch_len_interpolator.py:50-62`):

- Log-spaced near zero (5 points)
- Linear near zero (8 points)
- Quadratic from 0 to peak (n/3 points, dense near peak)
- Quadratic from peak to 3\*sigma (n/3 points)
- Quadratic from 3\*sigma to MAX_BRANCH_LENGTH (n/3 points, sparse tail)

This concentrates ~40 of 125 points near the peak while covering up to MAX_BRANCH_LENGTH = 4.0 subs/site.

## Fix (if GM tests fail at 1e-6)

Implement a non-uniform grid matching v0's multi-segment design, or use a two-tier approach: fine uniform grid over the peak region concatenated with a coarser grid for the tail. Requires `DistributionFunction` to support non-uniform grids or a resample step.

## Related issues

- Source: [kb/issues/M-timetree-branch-grid-uniform-resolution.md](../issues/M-timetree-branch-grid-uniform-resolution.md) -- delete after full resolution
