# Implement partition configuration layer

Four commands construct partitions inline with procedural code. The fix is not deduplicating init code but introducing a partition configuration layer so commands receive configured partitions.

## Current state

Each command (ancestral, optimize, timetree, prune) hardcodes partition creation: `create_fitch_partition` -> GTR resolution -> `log_gtr` -> `into_marginal`. Dense/sparse decided by `--dense` flag or `infer_dense()` stub. Always one sequence partition at index `0`. No `PartitionSpec` type exists.

## Target state

A `PartitionSpec` type describes what partition to create (alignment path, model, dense/sparse, partition index). A `create_partitions(specs, graph) -> Vec<Partition>` function builds partitions from specs. Configuration sources (CLI flags, config file, auto-detection) produce `Vec<PartitionSpec>`. Commands accept pre-built partitions.

## Implementation

1. Define `PartitionSpec { alignment, alphabet, model_source, dense_sparse, index }` in `partition/`
2. Implement `fn create_partitions(specs: &[PartitionSpec], graph: &Graph) -> Result<Vec<Partition>>` that handles create_fitch -> GTR resolution -> log -> convert
3. Update CLI arg structs to produce `Vec<PartitionSpec>` from `--aln`, `--model`, `--dense` flags (single-partition case maps to one-element vec)
4. Update all four commands to call `create_partitions` and receive results
5. Move command-specific post-init (write_gtr_json timing, rate normalization, trait-object wrapping) after the shared creation
6. This creates the extension point for future multi-partition sources (config file, `--segment`, auto-partitioning)

## Related issues

Source: [kb/issues/M-core-partition-init-orchestration-duplication.md](../issues/M-core-partition-init-orchestration-duplication.md)
Design: `kb/reports/partition-configuration.md`
Related: `kb/issues/H-core-command-module-shared-ops-entanglement.md`
Related: `kb/proposals/config-file-multi-partition.md`
