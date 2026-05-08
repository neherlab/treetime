# Promote pub(crate) optimize functions to pub for cross-crate access

## Description

Audit `pub(crate)` functions in `optimize/optimize_unified.rs` and promote the ones consumed by `commands/` to `pub`. Leave internal implementation details as `pub(crate)` or make them private.

## What to change

Promote to `pub`:

- `run_optimize_mixed_inner` - called by `commands/timetree/run.rs`
- `run_optimize_mixed_with_indel_rate` - called by optimize tests
- `evaluate_with_indels_log_lh_only` - called by `commands/timetree/inference/branch_length_likelihood.rs`
- `evaluate_with_indels` - called by Newton/Brent methods (internal to optimize/, but needs `pub` if methods move)

Keep `pub(crate)` or make private:

- `newton_tolerance_t`, `newton_tolerance_sqrt`, `newton_tolerance_log`
- `grid_search_inner`, `grid_search_branch_lengths`
- `reconcile_zero_boundary`, `is_zero_better_than_grid_best`
- `min_branch_length_for_indels`
- `GRID_SEARCH_MIN_UPPER`

## Related issues

- Source: [M-core-cli-library-separation-blockers](../issues/M-core-cli-library-separation-blockers.md) (B1)
