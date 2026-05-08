# Tighten pub(crate) visibility in optimize/

## Description

28 `pub(crate)` items in `optimize/`, most are internal helpers no code outside `optimize/` calls. Narrow to `pub(super)` or private. Promote true API functions to `pub`.

## Internal (narrow to pub(super) or private)

- `newton_tolerance_t`, `newton_tolerance_sqrt`, `newton_tolerance_log` (3 fns)
- `chain_rule_sqrt`, `chain_rule_log` (2 fns in method_newton.rs)
- `sqrt_step_lower_bound`, `log_step_lower_bound` (2 fns)
- `newton_inner`, `newton_sqrt_inner`, `newton_log_inner` (3 fns)
- `brent_inner`, `brent_sqrt_inner`, `brent_log_inner`, `brent_bracket` (4 fns)
- `NEWTON_REL_TOL`, `NEWTON_ABS_TOL`, `NEWTON_MAX_ITER`, `BRENT_MAX_ITER` (4 consts)
- `grid_search_inner`, `grid_search_branch_lengths` (2 fns)
- `reconcile_zero_boundary`, `is_zero_better_than_grid_best`, `min_branch_length_for_indels` (3 fns)
- `GRID_SEARCH_MIN_UPPER` (1 const)

## API (promote to pub)

- `run_optimize_mixed_inner` - called by `commands/timetree/run.rs`
- `run_optimize_mixed_with_indel_rate` - called by optimize tests
- `evaluate_with_indels_log_lh_only` - called by `commands/timetree/inference/branch_length_likelihood.rs`
- `evaluate_with_indels` - called by Newton/Brent methods

## Related issues

- Source: [M-core-remaining-architectural-debt-after-extraction](../issues/M-core-remaining-architectural-debt-after-extraction.md)
- Also tracked in: [M-core-cli-library-separation-blockers](../issues/M-core-cli-library-separation-blockers.md) (A7)
