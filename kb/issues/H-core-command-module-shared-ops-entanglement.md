# Command modules contain shared operations that belong in domain layers

Command modules (`commands/ancestral/`, `commands/clock/`, `commands/optimize/`, `commands/timetree/`, `commands/prune/`) accumulated domain logic that multiple commands and core layers depend on. This creates a reverse dependency: core modules (`representation/`, `gtr/`, `cli/`, `test_utils`) import from `commands/`, violating the expected layering where commands compose domain logic, not define it.

This initial sweep is not exhaustive. More investigation and more reorganization will be needed as the codebase evolves.

## Cross-command production imports

Five cross-command dependency edges exist in production code (excluding tests):

- ~~`clock` -> `ancestral`: `MethodAncestral`~~ **Resolved**: extracted to `commands/shared/args`
- ~~`clock` -> `timetree`: `BranchLengthMode`, `RerootMode`~~ **Resolved**: extracted to `commands/shared/args`
- ~~`timetree` -> `ancestral` (6 sites)~~ **Resolved**: ancestral algorithms extracted to top-level `src/ancestral/`; `get_common_length` moved to `src/seq/alignment`
- ~~`timetree` -> `clock` (20 sites)~~ **Resolved**: clock model domain extracted to top-level `src/clock/`
- ~~`timetree` -> `optimize` (8 sites)~~ **Resolved**: optimization machinery extracted to top-level `src/optimize/`; iteration helpers to `src/optimize/iteration.rs`
- ~~`optimize` -> `ancestral` (2 sites)~~ **Resolved**: ancestral algorithms extracted to top-level `src/ancestral/`; `get_common_length` moved to `src/seq/alignment`
- ~~`optimize` -> `prune` (1 site): `merge_shared_mutation_branches()`~~ **Resolved**: extracted to `partition/algo/topology_cleanup/merge_shared_mutations`
- `homoplasy` -> `ancestral`: `TreetimeAncestralArgs` - remaining cross-command arg dependency (minor: CLI struct embedding)

## Reverse dependencies (core importing from commands)

- ~~`partition/fitch_config.rs`, `partition/likelihood.rs` -> `get_common_length()` from `commands/ancestral/fitch`~~ **Resolved**: moved to `seq/alignment`
- ~~`PartitionOptimizeOps` from `commands/optimize/partition_ops`~~ **Resolved**: moved to `partition/traits.rs`
- ~~`PartitionRerootOps`, `PartitionTimetreeOps`, `PartitionTimetreeAll` from `commands/timetree/partition_ops`~~ **Resolved**: moved to `partition/traits.rs`
- ~~`RerootChanges` from `commands/clock/reroot`~~ **Resolved**: extracted to `partition/algo/topology_cleanup/reroot`
- ~~`ClockSet`, `ClockEdge`, `ClockNode`, `DateConstraintNode`, `TimetreeEdge`, `TimetreeNode` from `commands/clock/` and `commands/timetree/`~~ **Resolved**: `ClockSet` moved to `payload/clock_set.rs`; traits moved to `payload/traits.rs`
- ~~`OptimizationContribution` from `commands/optimize/optimize_unified`~~ **Resolved**: `OptimizationContribution` type + constructors moved to `partition/optimization_contribution.rs`; evaluation methods in `optimize/optimize_unified.rs` via split impl block
- ~~`ClockModel`, `ClockRegressionResult` from `commands/clock/`~~ **Resolved**: clock model domain extracted to top-level `src/clock/`
- ~~`compress_sequences()`, `get_common_length()`, `initialize_marginal()`, `update_marginal()` from `commands/ancestral/`~~ **Resolved**: ancestral algorithms extracted to top-level `src/ancestral/`; `get_common_length` moved to `src/seq/alignment`

## Resolved categories

### ~~Ancestral reconstruction~~ **Resolved**

Extracted to top-level `src/ancestral/` module: `fitch.rs`, `fitch_indel.rs`, `marginal.rs`. `get_common_length` moved to `src/seq/alignment.rs`. `commands/ancestral/` retains only `args.rs` and `run.rs`.

### ~~Clock model domain~~ **Resolved**

Extracted to top-level `src/clock/` module: `clock_model.rs`, `clock_regression.rs`, `clock_filter.rs`, `clock_output.rs`, `clock_graph.rs`, `date_constraints.rs`, `reroot.rs`, `rtt.rs`, `assign_dates.rs`, `find_best_root/`. `commands/clock/` retains only `args.rs` and `run.rs`.

### ~~Branch length optimization~~ **Resolved**

Extracted to top-level `src/optimize/` module: `optimize_unified.rs`, `optimize_indel.rs`, `method_brent.rs`, `method_newton.rs`, `optimize_eval.rs`, `optimize_dense_eval.rs`, `optimize_sparse_eval.rs`. Coefficient extraction types (`optimize_dense.rs`, `optimize_sparse.rs`, `optimization_contribution.rs`) moved to `partition/`. Iteration helpers (`save_branch_lengths`, `apply_damping`, topology cleanup) to `optimize/iteration.rs`. Domain arg types (`BranchOptMethod`, `InitialGuessMode`) to `optimize/args.rs`. `commands/optimize/` retains only `args.rs` and `run.rs`.

### ~~Partition and payload traits~~ **Resolved**

Traits and types moved to `representation/`: partition traits (`PartitionOptimizeOps`, `PartitionRerootOps`, `PartitionTimetreeOps`, `PartitionTimetreeAll`) to `partition/traits.rs`; payload traits (`ClockNode`, `ClockEdge`, `DateConstraintNode`, `TimetreeNode`, `TimetreeEdge`) to `payload/traits.rs`; `ClockSet` data type to `payload/clock_set.rs`.

### ~~Topology operations~~ **Resolved**

`merge_shared_mutation_branches()` extracted to `partition/algo/topology_cleanup/merge_shared_mutations`.

### ~~CLI args shared across commands~~ **Resolved**

`BranchLengthMode`, `RerootMode`, `MethodAncestral` extracted to `commands/shared/args.rs`. `BranchOptMethod`, `InitialGuessMode` extracted to `optimize/args.rs`.

## Remaining items

- `homoplasy` -> `ancestral`: `TreetimeAncestralArgs` CLI struct embedding (minor)
- Output helpers: `write_graph()` duplicated across `ancestral`, `optimize`, `prune` (tracked separately, not yet filed)
