# CLI code cannot be split into treetime-cli crate

After extracting domain modules (`ancestral/`, `clock/`, `optimize/`) from `commands/`, the `commands/` directory is thin enough to move into the `treetime-cli` crate. Three blockers prevent this.

## B1: (RESOLVED) pub(crate) functions consumed across the commands/ -> domain boundary

Resolved: `run_optimize_mixed_inner`, `run_optimize_mixed_with_indel_rate`, `evaluate_with_indels_log_lh_only`, `evaluate_with_indels` promoted to `pub` in `optimize/optimize_unified.rs`. Internal helpers (`newton_tolerance_*`, `grid_search_*`, `reconcile_zero_boundary`, `is_zero_better_than_grid_best`, `min_branch_length_for_indels`, `GRID_SEARCH_MIN_UPPER`) remain `pub(crate)`.

## B2: (RESOLVED) cli/ module serves commands only but lives in treetime crate

Resolved: `cli/rtt_chart.rs` and `cli/rtt_chart_render.rs` moved to `treetime-cli` crate. `run_clock` returns `ClockResult` and the CLI binary handles chart rendering. The `cli/` module and chart dependencies removed from `treetime` crate.

## B3: (RESOLVED) optimize/run.rs retains domain-ish functions consumed by domain tests

Resolved: `run_optimize_loop`, `collect_optimize_partitions`, `ConvergenceReason`, `OptimizeLoopResult` moved to `optimize/run_loop.rs`. `From<&BranchSplitArgs>` moved from `clock/find_best_root/params.rs` to `commands/clock/run.rs`.

## Non-blockers (already clean)

- `test_utils.rs` no longer imports from commands
- `cli/` no longer imports from commands
- `representation/` no longer imports from commands or optimize
- Domain modules (`ancestral/`, `optimize/`, `clock/`) are self-contained
