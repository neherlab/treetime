# CLI code cannot be split into treetime-cli crate

After extracting domain modules (`ancestral/`, `clock/`, `optimize/`) from `commands/`, the `commands/` directory is thin enough to move into the `treetime-cli` crate. Three blockers prevent this.

## B1: pub(crate) functions consumed across the commands/ -> domain boundary

`commands/timetree/run.rs` calls `run_optimize_mixed_inner` and `evaluate_with_indels_log_lh_only` from `optimize/optimize_unified.rs`, both `pub(crate)`. About 15 `pub(crate)` functions in `optimize/optimize_unified.rs` would need visibility audit:

Internal implementation details (should stay `pub(crate)` or become private):

- `newton_tolerance_t`, `newton_tolerance_sqrt`, `newton_tolerance_log`
- `grid_search_inner`, `grid_search_branch_lengths`
- `reconcile_zero_boundary`, `is_zero_better_than_grid_best`
- `min_branch_length_for_indels`
- `GRID_SEARCH_MIN_UPPER`

Genuine cross-module API (need `pub`):

- `run_optimize_mixed_inner` - called by `commands/timetree/run.rs`
- `run_optimize_mixed_with_indel_rate` - called by optimize tests
- `evaluate_with_indels_log_lh_only` - called by `commands/timetree/inference/branch_length_likelihood.rs`
- `evaluate_with_indels` - called by Newton/Brent methods

All iteration helpers in `optimize/iteration.rs` are already `pub` (moved during extraction).

## B2: cli/ module serves commands only but lives in treetime crate

`cli/rtt_chart.rs` and `cli/rtt_chart_render.rs` render root-to-tip regression charts. Called exclusively by `commands/clock/run.rs`. Import from `clock/` domain module (not from commands). If commands moves to `treetime-cli`, `cli/` goes with it or becomes public API.

Files: `packages/treetime/src/cli/rtt_chart.rs` (177 lines), `packages/treetime/src/cli/rtt_chart_render.rs` (174 lines).

## B3: optimize/run.rs retains domain-ish functions consumed by domain tests

`run_optimize_loop`, `collect_optimize_partitions`, `ConvergenceReason`, `OptimizeLoopResult` live in `commands/optimize/run.rs` but are consumed by tests in `optimize/__tests__/`. Moving commands to `treetime-cli` breaks these test imports.

Options:

- Move these functions to `optimize/` (they're optimization loop logic, arguably domain)
- Keep tests in `treetime-cli` alongside the functions (but they test domain behavior, not CLI)

Affected test files: `test_gm_optimize.rs`, `test_damping.rs`, `test_run_optimize_loop.rs`, `test_convergence_conditions.rs`, `test_convergence_sc2.rs`, `test_no_indels.rs`, `test_topology_cleanup.rs`, `test_dispatch_zero_boundary.rs`.

## Non-blockers (already clean)

- `test_utils.rs` no longer imports from commands
- `cli/` no longer imports from commands
- `representation/` no longer imports from commands or optimize
- Domain modules (`ancestral/`, `optimize/`, `clock/`) are self-contained
