# Branch distribution grid uses uniform spacing

The branch distribution grid in `create_simple_grid()` uses `Array1::linspace` (uniform spacing) over `[one_mutation * 0.01, min(center * 5, MAX_BRANCH_LENGTH)]`. The grid is now branch-length-informed (extent scales with the ML branch length, capped at 5.0 subs/site) and uses 300 points, a substantial improvement over the prior 1000-point grid that spread resolution across a 200-year time window.

## Current status

The branch-length-informed grid concentrates points near the peak for typical branches. The resolution improvement has not been empirically validated: the golden master tests (`test_gm_runner_marginal_dense`, `test_gm_runner_marginal_sparse`, `test_gm_runner_poisson`) remain `#[ignore]`d at target tolerance `epsilon = 1e-6`. Un-ignoring them is the acceptance criterion for closing this issue.

## Remaining concern

For very short branches where `one_mutation * 10` governs the extent, the uniform grid may still under-resolve the peak relative to v0's non-uniform grid. A non-uniform grid concentrating points near the peak while maintaining broad tail coverage (as v0 does) would further improve resolution for these cases.

## v0 reference

v0 uses a 5-segment non-uniform grid (`branch_len_interpolator.py:50-62`):

- Log-spaced near zero (5 points)
- Linear near zero (8 points)
- Quadratic from 0 to peak (n/3 points, dense near peak)
- Quadratic from peak to 3\*sigma (n/3 points)
- Quadratic from 3\*sigma to MAX_BRANCH_LENGTH (n/3 points, sparse tail)

This concentrates ~40 of 125 points near the peak while covering up to MAX_BRANCH_LENGTH = 4.0 subs/site.

## Related tickets

- [kb/tickets/timetree-switch-branch-grid-to-nonuniform-spacing.md](../tickets/timetree-switch-branch-grid-to-nonuniform-spacing.md)
