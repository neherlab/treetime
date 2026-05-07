# Extract optimization machinery from commands/optimize/ to domain module

## Description

`commands/optimize/` contains shared branch length optimization machinery (5 files, 1433 lines) consumed by `timetree` (8 import sites) and `representation/partition/`. These are domain algorithms, not command logic.

## What to move

5 files in `commands/optimize/`:

- `optimize_unified.rs` (806 lines) -- `run_optimize_mixed()`, `OptimizationContribution`, `evaluate_with_indels_log_lh_only()`
- `optimize_indel.rs` (188 lines) -- `estimate_indel_rate()`
- `method_brent.rs` (185 lines) -- Brent optimization method
- `method_newton.rs` (233 lines) -- Newton optimization method
- `partition_ops.rs` (21 lines) -- `PartitionOptimizeOps`

Total: 1433 lines.

## Consumers

### timetree (8 import sites)

`run_optimize_mixed()`, `OptimizationContribution`, `evaluate_with_indels_log_lh_only()`, `estimate_indel_rate()`, `BranchOptMethod`, `PartitionOptimizeOps`, `apply_damping()`, `save_branch_lengths()`

### representation/partition/ (reverse dependency -- core importing from commands)

- `representation/partition/marginal_dense.rs#L2` -> `PartitionOptimizeOps` from `commands/optimize/partition_ops`
- `representation/partition/marginal_sparse.rs#L3` -> `PartitionOptimizeOps` from `commands/optimize/partition_ops`

## Target location

Extract to a new domain module outside `commands/` (e.g. `src/optimize/` or a new crate). The `commands/optimize/` command module should become a thin wrapper.

## Related issues

- Source: [H-core-command-module-shared-ops-entanglement.md](../issues/H-core-command-module-shared-ops-entanglement.md) -- delete after full resolution
