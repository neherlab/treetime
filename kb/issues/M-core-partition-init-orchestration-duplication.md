# Partition creation is hardcoded per-command instead of configured

## Current state

`partition/create.rs` provides `create_marginal_partition()` consolidating the 3-way branch (sparse, dense+infer, dense+named) without `write_gtr_json` during init. All pipeline functions use it. The full `PartitionSpec`/config layer for multi-partition support is not yet implemented.

## Original description

Four commands construct partitions inline with procedural code: `create_fitch_partition` -> GTR resolution -> `log_gtr` -> `into_marginal_sparse`/`into_marginal_dense`. Each command implements its own variant with different models, different dense/sparse decisions, different `write_gtr_json` placement, and different post-init handling. These are not duplicated identical partitions but independent partition setups that share no configuration layer.

The multi-partition infrastructure (`Vec<Arc<RwLock<dyn PartitionTimetreeAll>>>`) exists and works (mugration adds a discrete partition), but no command creates more than one sequence partition. No `PartitionConfig` or `PartitionSpec` type exists. Partition creation is entangled with command orchestration.

## Locations

- `packages/treetime/src/commands/ancestral/run.rs:134-240:` three paths (sparse, dense+infer, dense+named) -- `write_gtr_json` after reconstruction (sparse/dense+infer) or before (dense+named)
- `packages/treetime/src/commands/optimize/run.rs:67-117:` three paths -- `write_gtr_json` deferred until after `normalize_partition_rates`
- `packages/treetime/src/commands/timetree/initialization.rs:86-124:` three paths -- `write_gtr_json` immediately after partition creation, wraps into `Arc<RwLock<dyn PartitionTimetreeAll>>`
- `packages/treetime/src/commands/prune/run.rs:48-60:` sparse-only path -- `write_gtr_json` before `into_marginal_sparse`, hardcoded JC69 model

## Impact

Adding multi-partition support (multi-segment genomes, codon-position partitioning, config file) requires changing all four commands individually. The `write_gtr_json` placement inconsistency means GTR JSON output reflects different states per command.

## Action

Implement partition configuration: commands receive configured partitions rather than constructing them inline. See `kb/reports/partition-configuration.md` for the full design landscape covering CLI flags, config files, auto-detection heuristics, and model binding.

## Related

- `kb/issues/H-core-command-module-shared-ops-entanglement.md` -- broader command module entanglement
- `kb/issues/N-optimize-multi-alignment-input.md` -- optimize accepts only single alignment
- `kb/issues/N-io-multi-segment-genome-input.md` -- multi-segment genome input not wired
- `kb/issues/N-representation-infer-dense-stub.md` -- dense/sparse auto-detection stub
- `kb/proposals/config-file-multi-partition.md` -- config file proposal
- `kb/reports/partition-configuration.md` -- full partition configuration report

## Related tickets

- [kb/tickets/architecture-extract-shared-partition-init-helper.md](../tickets/architecture-extract-shared-partition-init-helper.md)
