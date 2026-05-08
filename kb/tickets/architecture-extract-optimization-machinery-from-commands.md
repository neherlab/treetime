# Extract optimization machinery from commands/optimize/ to domain module

## Description

`commands/optimize/` contains shared branch length optimization machinery (5 files, 1433 lines) consumed by `timetree` (8 import sites) and `representation/partition/`. These are domain algorithms, not command logic.

## What to move

4 files in `commands/optimize/`:

- `optimize_unified.rs` (806 lines) -- `run_optimize_mixed()`, `OptimizationContribution`, `evaluate_with_indels_log_lh_only()`
- `optimize_indel.rs` (188 lines) -- `estimate_indel_rate()`
- `method_brent.rs` (185 lines) -- Brent optimization method
- `method_newton.rs` (233 lines) -- Newton optimization method

Total: 1412 lines.

`PartitionOptimizeOps` already moved to `representation/partition/traits.rs`.

## Consumers

### timetree (8 import sites)

`run_optimize_mixed()`, `OptimizationContribution`, `evaluate_with_indels_log_lh_only()`, `estimate_indel_rate()`, `BranchOptMethod`, `PartitionOptimizeOps`, `apply_damping()`, `save_branch_lengths()`

### representation/partition/ (residual reverse dependency)

- `representation/partition/marginal_dense.rs` -> `OptimizationContribution` from `commands/optimize/optimize_unified`
- `representation/partition/marginal_sparse.rs` -> `OptimizationContribution` from `commands/optimize/optimize_unified`
- `representation/partition/traits.rs` -> `OptimizationContribution` from `commands/optimize/optimize_unified` (return type of `PartitionOptimizeOps::create_edge_contribution`)

## Target location

Extract to a new domain module outside `commands/` (e.g. `src/optimize/` or a new crate). The `commands/optimize/` command module should become a thin wrapper.

## Related issues

- Source: [H-core-command-module-shared-ops-entanglement.md](../issues/H-core-command-module-shared-ops-entanglement.md) -- delete after full resolution
