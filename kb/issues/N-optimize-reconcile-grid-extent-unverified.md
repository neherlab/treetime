# `reconcile_zero_boundary` grid extent argument is unverified

## Problem

`reconcile_zero_boundary` at [packages/treetime/src/optimize/dispatch.rs](../../packages/treetime/src/optimize/dispatch.rs) takes a `branch_length_for_extent` parameter separate from `candidate`. The separation exists because `candidate` may be clamped to $0$ or to a tiny floor like $10^{-12}$ after the inner solver runs, and the grid extent must remain tied to a meaningful scale. The dispatcher in `run_optimize_mixed` passes the edge's INPUT branch length (the value before optimization) as the extent.

This choice is correct in the common case: the input branch length is the optimizer's best prior estimate of the true scale. The corner case where the input is uninformative or wildly misleading is not covered:

- Input branch length is zero (e.g. uncalibrated tree): grid extent falls back to `max(0.1 * one_mutation, GRID_SEARCH_MIN_UPPER = 0.5)`, which may not contain the true mode if the true scale is outside $[0.1 \cdot \text{one\_mutation}, 0.5]$.
- Input branch length is much smaller than the true mode (e.g. bad initial guess): `1.5 * bl + one_mutation < true_mode`, and the cap at `0.5` saves the day only if `true_mode < 0.5`.
- Input branch length is much larger than the true mode: grid scans a range far from the true mode, may miss it entirely if the grid's log spacing puts all points above the true scale.

## Impact

Negligible when the initial guess is calibrated (the standard workflow). On datasets with uncalibrated input trees or bad initial guesses combined with non-unimodal models that force `reconcile_zero_boundary` to grid-search, the reconciled branch length lands on whichever grid point beats zero, and there is no test that proves the true mode is on the grid. The behavior is unverified.

## Proposed action

1. Test with a fixture where the input branch length is $100\times$ away from the true mode in each direction (above and below), verify the grid still covers the true mode.
2. If the test fails, consider passing a more informative extent: e.g. the maximum of the input branch length and `one_mutation * sqrt(sequence_length)` (the standard deviation of the Hamming distance distribution under a neutral model), or an invariant derived from the current contribution coefficients.
3. Document the admissibility contract: "the caller must pass a `branch_length_for_extent` whose true optimum lies within $[\text{extent} / 10, 1.5 \cdot \text{extent} + \text{one\_mutation}]$ or $[\text{extent} / 10, 0.5]$, whichever is larger".

## Cross-references

- [packages/treetime/src/optimize/dispatch.rs](../../packages/treetime/src/optimize/dispatch.rs) (`reconcile_zero_boundary`, `grid_search_branch_lengths`)
