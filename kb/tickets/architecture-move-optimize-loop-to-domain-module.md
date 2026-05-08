# Move run_optimize_loop and related functions from commands/ to optimize/

## Description

`run_optimize_loop`, `collect_optimize_partitions`, `ConvergenceReason`, `OptimizeLoopResult` live in `commands/optimize/run.rs` but are consumed by 8 test files in `optimize/__tests__/`. These are optimization loop logic (convergence detection, partition collection), not CLI orchestration. Move them to `optimize/` so tests and the functions colocate in the same crate.

## What to change

Move from `commands/optimize/run.rs` to `optimize/iteration.rs` (or a new `optimize/loop.rs`):

- `fn run_optimize_loop` (the iteration loop with convergence detection)
- `fn collect_optimize_partitions` (partition type erasure helper)
- `fn compute_iteration_likelihood` (private helper called by run_optimize_loop)
- `enum ConvergenceReason`
- `struct OptimizeLoopResult`

After move, `commands/optimize/run.rs` retains only:

- `fn run_optimize` (CLI entry point: parse args, read files, call domain functions, write output)
- `struct TreetimeOptimizeParams`

## Test files that currently cross the boundary

- `test_gm_optimize.rs`
- `test_damping.rs`
- `test_run_optimize_loop.rs`
- `test_convergence_conditions.rs`
- `test_convergence_sc2.rs`
- `test_no_indels.rs`
- `test_topology_cleanup.rs`
- `test_dispatch_zero_boundary.rs`

## Related issues

- Source: [M-core-cli-library-separation-blockers](../issues/M-core-cli-library-separation-blockers.md) (B3)
