# Command modules contain shared operations that belong in domain layers

Command modules (`commands/ancestral/`, `commands/clock/`, `commands/optimize/`, `commands/timetree/`, `commands/prune/`) accumulated domain logic that multiple commands and core layers depend on. This creates a reverse dependency: core modules (`representation/`, `gtr/`, `cli/`, `test_utils`) import from `commands/`, violating the expected layering where commands compose domain logic, not define it.

This initial sweep is not exhaustive. More investigation and more reorganization will be needed as the codebase evolves.

## Cross-command production imports

Six cross-command dependency edges exist in production code (excluding tests):

- `clock` -> `ancestral`: `MethodAncestral` (`#MethodAncestral`) [packages/treetime/src/commands/ancestral/args.rs#L11](../../packages/treetime/src/commands/ancestral/args.rs#L11)
- `clock` -> `timetree`: `BranchLengthMode` (`#BranchLengthMode`) [packages/treetime/src/commands/timetree/args.rs#L21](../../packages/treetime/src/commands/timetree/args.rs#L21), `RerootMode` (`#RerootMode`) [packages/treetime/src/commands/timetree/args.rs#L40](../../packages/treetime/src/commands/timetree/args.rs#L40)
- `timetree` -> `ancestral` (6 sites): `compress_sequences()` (`#compress_sequences`) [packages/treetime/src/commands/ancestral/fitch.rs#L526](../../packages/treetime/src/commands/ancestral/fitch.rs#L526), `get_common_length()` (`#get_common_length`) [packages/treetime/src/commands/ancestral/fitch.rs#L616](../../packages/treetime/src/commands/ancestral/fitch.rs#L616), `initialize_marginal()` (`#initialize_marginal`) [packages/treetime/src/commands/ancestral/marginal.rs#L21](../../packages/treetime/src/commands/ancestral/marginal.rs#L21), `update_marginal()` (`#update_marginal`) [packages/treetime/src/commands/ancestral/marginal.rs#L41](../../packages/treetime/src/commands/ancestral/marginal.rs#L41)
- `timetree` -> `clock` (20 sites): `ClockModel` (`#ClockModel`) [packages/treetime/src/commands/clock/clock_model.rs#L31](../../packages/treetime/src/commands/clock/clock_model.rs#L31), `ClockSet` (`#ClockSet`) [packages/treetime/src/commands/clock/clock_set.rs#L9](../../packages/treetime/src/commands/clock/clock_set.rs#L9), `ClockNode` (`#ClockNode`) [packages/treetime/src/commands/clock/clock_traits.rs#L5](../../packages/treetime/src/commands/clock/clock_traits.rs#L5), `ClockEdge` (`#ClockEdge`) [packages/treetime/src/commands/clock/clock_traits.rs#L16](../../packages/treetime/src/commands/clock/clock_traits.rs#L16), `ClockParams` (`#ClockParams`) [packages/treetime/src/commands/clock/clock_regression.rs#L19](../../packages/treetime/src/commands/clock/clock_regression.rs#L19), `RerootChanges` (`#RerootChanges`) [packages/treetime/src/commands/clock/reroot.rs#L78](../../packages/treetime/src/commands/clock/reroot.rs#L78), `RerootParams` (`#RerootParams`) [packages/treetime/src/commands/clock/reroot.rs#L17](../../packages/treetime/src/commands/clock/reroot.rs#L17), `BranchPointOptimizationParams` (`#BranchPointOptimizationParams`), `clock_filter_inplace()` (`#clock_filter_inplace`) [packages/treetime/src/commands/clock/clock_filter.rs#L21](../../packages/treetime/src/commands/clock/clock_filter.rs#L21), `estimate_clock_model_with_reroot()` (`#estimate_clock_model_with_reroot`) [packages/treetime/src/commands/clock/clock_regression.rs#L131](../../packages/treetime/src/commands/clock/clock_regression.rs#L131), `estimate_clock_model_with_reroot_policy()` (`#estimate_clock_model_with_reroot_policy`) [packages/treetime/src/commands/clock/clock_regression.rs#L161](../../packages/treetime/src/commands/clock/clock_regression.rs#L161), `load_date_constraints()` (`#load_date_constraints`) [packages/treetime/src/commands/clock/date_constraints.rs#L23](../../packages/treetime/src/commands/clock/date_constraints.rs#L23), `write_clock_model()` (`#write_clock_model`) [packages/treetime/src/commands/clock/clock_output.rs#L11](../../packages/treetime/src/commands/clock/clock_output.rs#L11)
- `timetree` -> `optimize` (8 sites): `run_optimize_mixed()` (`#run_optimize_mixed`) [packages/treetime/src/commands/optimize/optimize_unified.rs#L533](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L533), `OptimizationContribution` (`#OptimizationContribution`) [packages/treetime/src/commands/optimize/optimize_unified.rs#L96](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L96), `evaluate_with_indels_log_lh_only()` (`#evaluate_with_indels_log_lh_only`) [packages/treetime/src/commands/optimize/optimize_unified.rs#L259](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L259), `estimate_indel_rate()` (`#estimate_indel_rate`) [packages/treetime/src/commands/optimize/optimize_indel.rs#L55](../../packages/treetime/src/commands/optimize/optimize_indel.rs#L55), `BranchOptMethod` (`#BranchOptMethod`), `PartitionOptimizeOps` (`#PartitionOptimizeOps`) [packages/treetime/src/commands/optimize/partition_ops.rs#L13](../../packages/treetime/src/commands/optimize/partition_ops.rs#L13), `apply_damping()` (`#apply_damping`) [packages/treetime/src/commands/optimize/run.rs#L513](../../packages/treetime/src/commands/optimize/run.rs#L513), `save_branch_lengths()` (`#save_branch_lengths`) [packages/treetime/src/commands/optimize/run.rs#L460](../../packages/treetime/src/commands/optimize/run.rs#L460)
- `optimize` -> `ancestral` (2 sites): `compress_sequences()`, `get_common_length()`, `initialize_marginal()`, `update_marginal()` (same definitions as above)
- `optimize` -> `prune` (1 site): `merge_shared_mutation_branches()` (`#merge_shared_mutation_branches`) [packages/treetime/src/commands/prune/run.rs#L294](../../packages/treetime/src/commands/prune/run.rs#L294)
- `homoplasy` -> `ancestral`: `TreetimeAncestralArgs` (`#TreetimeAncestralArgs`) [packages/treetime/src/commands/ancestral/args.rs#L19](../../packages/treetime/src/commands/ancestral/args.rs#L19)

## Reverse dependencies (core importing from commands)

Core modules import from `commands/`, inverting the expected dependency direction:

- [packages/treetime/src/representation/partition/fitch_config.rs#L2](../../packages/treetime/src/representation/partition/fitch_config.rs#L2), [packages/treetime/src/representation/partition/likelihood.rs#L2](../../packages/treetime/src/representation/partition/likelihood.rs#L2) -> `get_common_length()` from `commands/ancestral/fitch`
- [packages/treetime/src/representation/partition/marginal_dense.rs#L2](../../packages/treetime/src/representation/partition/marginal_dense.rs#L2), [packages/treetime/src/representation/partition/marginal_sparse.rs#L3](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L3) -> `PartitionOptimizeOps` from `commands/optimize/partition_ops`
- [packages/treetime/src/representation/partition/marginal_dense.rs#L3](../../packages/treetime/src/representation/partition/marginal_dense.rs#L3), [packages/treetime/src/representation/partition/marginal_sparse.rs#L4](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L4), [packages/treetime/src/representation/partition/timetree.rs#L1](../../packages/treetime/src/representation/partition/timetree.rs#L1) -> `PartitionRerootOps`, `PartitionTimetreeAll` from `commands/timetree/partition_ops`
- [packages/treetime/src/representation/partition/marginal_sparse.rs#L2](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L2) -> `RerootChanges` from `commands/clock/reroot`
- [packages/treetime/src/representation/payload/timetree.rs#L1-L4](../../packages/treetime/src/representation/payload/timetree.rs#L1-L4) -> `ClockSet`, `ClockEdge`, `ClockNode`, `DateConstraintNode`, `TimetreeEdge`, `TimetreeNode` from `commands/clock/` and `commands/timetree/`
- [packages/treetime/src/cli/rtt_chart.rs#L2-L3](../../packages/treetime/src/cli/rtt_chart.rs#L2-L3), [packages/treetime/src/cli/rtt_chart_render.rs#L1-L2](../../packages/treetime/src/cli/rtt_chart_render.rs#L1-L2) -> `ClockModel`, `ClockRegressionResult` from `commands/clock/`
- [packages/treetime/src/test_utils.rs#L2-L3](../../packages/treetime/src/test_utils.rs#L2-L3) -> `compress_sequences()`, `get_common_length()`, `initialize_marginal()`, `update_marginal()` from `commands/ancestral/`

## Affected modules by category

### Ancestral reconstruction (797 lines)

[packages/treetime/src/commands/ancestral/fitch.rs](../../packages/treetime/src/commands/ancestral/fitch.rs) (640 lines) and [packages/treetime/src/commands/ancestral/marginal.rs](../../packages/treetime/src/commands/ancestral/marginal.rs) (157 lines) define the core reconstruction algorithms. Consumed by `ancestral`, `optimize`, `timetree`, `prune`, `representation/`, `gtr/`, and `test_utils`. These are domain algorithms, not command logic.

### Clock model (1443 lines)

`commands/clock/` contains 11 non-command files: [packages/treetime/src/commands/clock/clock_model.rs](../../packages/treetime/src/commands/clock/clock_model.rs), [packages/treetime/src/commands/clock/clock_set.rs](../../packages/treetime/src/commands/clock/clock_set.rs), [packages/treetime/src/commands/clock/clock_traits.rs](../../packages/treetime/src/commands/clock/clock_traits.rs), [packages/treetime/src/commands/clock/clock_regression.rs](../../packages/treetime/src/commands/clock/clock_regression.rs), [packages/treetime/src/commands/clock/clock_filter.rs](../../packages/treetime/src/commands/clock/clock_filter.rs), [packages/treetime/src/commands/clock/clock_output.rs](../../packages/treetime/src/commands/clock/clock_output.rs), [packages/treetime/src/commands/clock/reroot.rs](../../packages/treetime/src/commands/clock/reroot.rs), [packages/treetime/src/commands/clock/rtt.rs](../../packages/treetime/src/commands/clock/rtt.rs), [packages/treetime/src/commands/clock/date_constraints.rs](../../packages/treetime/src/commands/clock/date_constraints.rs), [packages/treetime/src/commands/clock/clock_graph.rs](../../packages/treetime/src/commands/clock/clock_graph.rs), [packages/treetime/src/commands/clock/assign_dates.rs](../../packages/treetime/src/commands/clock/assign_dates.rs). Consumed by `timetree` (20 import sites), `representation/payload/`, `representation/partition/`, and `cli/`. These define the clock model domain, not the `clock` command.

### Branch length optimization (1433 lines)

`commands/optimize/` contains shared optimization machinery: [packages/treetime/src/commands/optimize/optimize_unified.rs](../../packages/treetime/src/commands/optimize/optimize_unified.rs) (806 lines), [packages/treetime/src/commands/optimize/optimize_indel.rs](../../packages/treetime/src/commands/optimize/optimize_indel.rs) (188 lines), [packages/treetime/src/commands/optimize/method_brent.rs](../../packages/treetime/src/commands/optimize/method_brent.rs) (185 lines), [packages/treetime/src/commands/optimize/method_newton.rs](../../packages/treetime/src/commands/optimize/method_newton.rs) (233 lines), [packages/treetime/src/commands/optimize/partition_ops.rs](../../packages/treetime/src/commands/optimize/partition_ops.rs) (21 lines). Consumed by `timetree` (8 import sites) and `representation/partition/`. The optimization algorithms are domain logic used across commands.

### Partition and payload traits (78 lines)

[packages/treetime/src/commands/timetree/partition_ops.rs](../../packages/treetime/src/commands/timetree/partition_ops.rs) (60 lines) and [packages/treetime/src/commands/timetree/timetree_traits.rs](../../packages/treetime/src/commands/timetree/timetree_traits.rs) (18 lines) define traits consumed by `representation/partition/` and `representation/payload/`. These are trait definitions for the representation layer, placed in a command module.

### Topology operations

`merge_shared_mutation_branches()` (`#merge_shared_mutation_branches`) in [packages/treetime/src/commands/prune/run.rs#L294](../../packages/treetime/src/commands/prune/run.rs#L294) is consumed by `commands/optimize/`. Already tracked in [L-topology-cleanup-move-merge-shared-mutations](L-topology-cleanup-move-merge-shared-mutations.md).

### CLI args shared across commands

[packages/treetime/src/commands/timetree/args.rs](../../packages/treetime/src/commands/timetree/args.rs) exports `BranchLengthMode` (`#BranchLengthMode`, L21) and `RerootMode` (`#RerootMode`, L40) consumed by [packages/treetime/src/commands/clock/args.rs](../../packages/treetime/src/commands/clock/args.rs). [packages/treetime/src/commands/ancestral/args.rs](../../packages/treetime/src/commands/ancestral/args.rs) exports `MethodAncestral` (`#MethodAncestral`, L11) consumed by [packages/treetime/src/commands/clock/args.rs](../../packages/treetime/src/commands/clock/args.rs) and [packages/treetime/src/commands/timetree/args.rs](../../packages/treetime/src/commands/timetree/args.rs).

### Output helpers

`write_graph()` duplicated across `ancestral`, `optimize`, `prune`. Already tracked in [L-core-duplicate-write-graph](L-core-duplicate-write-graph.md).

## Proposed reorganization

The command modules should be thin wrappers: parse CLI args, call domain functions, write output. Domain logic should live in domain layers.

Candidate extraction targets (locations are proposals, not final):

1. **Ancestral reconstruction algorithms** ([packages/treetime/src/commands/ancestral/fitch.rs](../../packages/treetime/src/commands/ancestral/fitch.rs), [packages/treetime/src/commands/ancestral/marginal.rs](../../packages/treetime/src/commands/ancestral/marginal.rs)) -> a new domain module outside `commands/` (e.g. `src/ancestral/` or a new crate)
2. **Clock model domain** (11 files, 1443 lines in `commands/clock/`) -> a new domain module outside `commands/` (e.g. `src/clock/` or a new crate)
3. **Branch length optimization** ([packages/treetime/src/commands/optimize/optimize_unified.rs](../../packages/treetime/src/commands/optimize/optimize_unified.rs), [packages/treetime/src/commands/optimize/optimize_indel.rs](../../packages/treetime/src/commands/optimize/optimize_indel.rs), [packages/treetime/src/commands/optimize/method_brent.rs](../../packages/treetime/src/commands/optimize/method_brent.rs), [packages/treetime/src/commands/optimize/method_newton.rs](../../packages/treetime/src/commands/optimize/method_newton.rs), [packages/treetime/src/commands/optimize/partition_ops.rs](../../packages/treetime/src/commands/optimize/partition_ops.rs)) -> a new domain module outside `commands/` (e.g. `src/optimize/` or a new crate)
4. **Partition/payload traits** ([packages/treetime/src/commands/timetree/partition_ops.rs](../../packages/treetime/src/commands/timetree/partition_ops.rs), [packages/treetime/src/commands/timetree/timetree_traits.rs](../../packages/treetime/src/commands/timetree/timetree_traits.rs)) -> `representation/` where the implementing types live
5. **Shared CLI arg types** (`BranchLengthMode`, `RerootMode`, `MethodAncestral`) -> a shared args module or `treetime-cli`
6. **I/O helpers** (`write_graph()`, `write_clock_model()`) -> `treetime-io`

Each extraction should be done incrementally. The dependency graph is DAG-shaped (no cycles between domain clusters), so extractions can proceed independently.

## Related

- [L-core-duplicate-write-graph](L-core-duplicate-write-graph.md) - output helper duplication (subset of this issue)
- [L-topology-cleanup-move-merge-shared-mutations](L-topology-cleanup-move-merge-shared-mutations.md) - topology operation misplacement (subset of this issue)
