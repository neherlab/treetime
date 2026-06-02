# Make the partition the single source of truth for the pipeline GTR

## Task

Eliminate the duplicate GTR carried by `PartitionCreated.gtr` so every pipeline derives its output GTR from the owning partition at output time, after all in-place mutations (normalization, refinement).

## Steps

1. optimize (`packages/treetime/src/optimize/pipeline.rs`): stop using `created.gtr`. After the optimize loop and final `update_marginal`, read the GTR from the live partition (the non-empty sparse or dense partition) via `HasGtr::gtr()` and return it as `OptimizeOutput.gtr`. This restores the normalized `mu` (== 1.0 for a single partition) and matches v0.
2. timetree (`packages/treetime/src/timetree/pipeline.rs`, `initialize_partitions_from_params`): read the GTR from `created.partition` (the owning `MarginalPartition`) instead of `created.gtr`. Behavior-preserving today, since timetree does not mutate the GTR.
3. `partition/create.rs`: remove the `gtr` field from `PartitionCreated` (keep the local `gtr` binding for `log_gtr` and `into_marginal_*`). The bug becomes unrepresentable: there is no standalone GTR left to go stale.

## Test

Add a unit test asserting `optimize --gtr=infer` produces an `OptimizeOutput.gtr` with `mu` normalized to 1.0 (the post-condition of `normalize_partition_rates` for a single partition). The raw inferred `mu` is essentially never 1.0, so the test fails on the stale snapshot ("see red"). Use a small real dataset (`data/flu/h3n2/20`).

## Out of scope

- Removing unwired-but-intentional scaffolding (`TimetreeParams` fields for not-yet-parted features, auspice `AuspiceWrite::new` trait stub).
- A `MarginalPartition::gtr()` accessor or `PartitionTimetreeAll` trait extension (not required; callers read the concrete partition type).

## Related issues

- Source: [H-core-pipeline-gtr-stale-snapshot.md](../issues/H-core-pipeline-gtr-stale-snapshot.md) -- delete after full resolution
- [M-core-partition-init-orchestration-duplication.md](../issues/M-core-partition-init-orchestration-duplication.md) -- same `create_marginal_partition` area
