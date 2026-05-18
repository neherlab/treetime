# optimize/ internal file cohesion

Two files in `optimize/` combine multiple responsibilities that would be clearer as separate modules.

## optimize_unified.rs (701 lines, 12 public/pub(crate) items)

Three distinct responsibilities in one file:

1. Data types and evaluation combinators: `OptimizationMetrics`, `evaluate_mixed`, `evaluate_with_indels`, `evaluate_mixed_log_lh_only`, `evaluate_with_indels_log_lh_only`
2. Zero-boundary decision logic: `is_zero_branch_optimal`, `is_zero_better_than_grid_best`, `reconcile_zero_boundary`, `min_branch_length_for_indels` (~130 lines of dense doc comments explaining multimodal surface edge cases)
3. Per-edge dispatch and initial guess: `run_optimize_mixed`, `run_optimize_mixed_inner`, `run_optimize_mixed_with_indel_rate`, `initial_guess_mixed`

Grid search (`grid_search_branch_lengths`, `grid_search_inner`) and Newton tolerance functions (`newton_tolerance_t`, `newton_tolerance_sqrt`, `newton_tolerance_log`) are also here, loosely coupled to group 3.

## iteration.rs (336 lines, 10 public items)

Mixed abstraction levels under an opaque name:

- Loop bookkeeping: `save_branch_lengths`, `restore_branch_lengths`, `apply_damping` (branch-length save/restore/damping between iterations)
- Topology surgery: `prune_and_merge_in_loop`, `find_zero_optimal_internal_edges` (calls into `representation/algo/topology_cleanup/`)
- Setup/initialization: `apply_initial_guess_mode`, `normalize_partition_rates`, `any_edge_missing_branch_length`, `any_indel_edge_has_zero_branch_length`

"Iteration" describes when these run (during the loop), not what they do.

## Naming inconsistencies

- `optimize/` vs `timetree/optimization/` for sibling modules
- `optimize_eval.rs`, `optimize_unified.rs`, `optimize_indel.rs` stutter the parent module name. `method_brent.rs` and `method_newton.rs` use a different prefix convention. No consistent scheme within the module
- `args.rs` contains domain enums (`BranchOptMethod`, `InitialGuessMode`) used deep in the optimizer, not CLI argument types

## Impact

Negligible. No correctness or performance issue. Navigation friction when reading or modifying the optimizer.
