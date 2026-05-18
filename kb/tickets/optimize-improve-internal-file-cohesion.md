# Improve optimize/ internal file cohesion

## Description

Split `optimize_unified.rs` and `iteration.rs` along responsibility boundaries. Rename files to remove stutter and apply consistent naming.

## Proposed splits

### optimize_unified.rs -> 3 files

1. `likelihood.rs` - `OptimizationMetrics` struct, `evaluate_mixed`, `evaluate_with_indels`, log-likelihood-only variants
2. `zero_boundary.rs` - `is_zero_branch_optimal`, `is_zero_better_than_grid_best`, `reconcile_zero_boundary`, `min_branch_length_for_indels`, grid search functions
3. `dispatch.rs` - `run_optimize_mixed`, `run_optimize_mixed_inner`, `run_optimize_mixed_with_indel_rate`, `initial_guess_mixed`

Newton tolerance functions (`newton_tolerance_t`, `newton_tolerance_sqrt`, `newton_tolerance_log`) move to `method_newton.rs` alongside the constants they reference.

### iteration.rs cleanup

1. Keep `iteration.rs` with branch-length bookkeeping: `save_branch_lengths`, `restore_branch_lengths`, `apply_damping`, `DAMPING_FLOOR`
2. Move `prune_and_merge_in_loop`, `find_zero_optimal_internal_edges` to `run_loop.rs` (only caller)
3. Move `apply_initial_guess_mode`, `normalize_partition_rates`, `any_edge_missing_branch_length`, `any_indel_edge_has_zero_branch_length` to `run_loop.rs` (only production caller)

### Naming cleanup

- Remove `optimize_` prefix from `optimize_eval.rs`, `optimize_indel.rs`, `optimize_dense_eval.rs`, `optimize_sparse_eval.rs`
- Rename `args.rs` to `params.rs` (contents are domain enums, not CLI args)

## Constraints

- Pure refactoring: no behavior change, no public API change beyond module paths
- All `pub(crate)` visibility stays the same (narrowing tracked separately in `architecture-tighten-optimize-visibility`)
- Test files stay in `__tests__/` - update import paths only

## Related issues

- Source: [N-optimize-internal-file-cohesion](../issues/N-optimize-internal-file-cohesion.md)
- Adjacent: [architecture-tighten-optimize-visibility](architecture-tighten-optimize-visibility.md) (visibility narrowing, orthogonal)
- Adjacent: [M-core-remaining-architectural-debt-after-extraction](../issues/M-core-remaining-architectural-debt-after-extraction.md) (inter-module concerns, separate scope)
