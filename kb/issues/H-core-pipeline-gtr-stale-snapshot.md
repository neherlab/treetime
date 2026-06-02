# Pipeline serializes a stale GTR snapshot instead of the partition's live model

The shared partition constructor `create_marginal_partition` returns a standalone clone of the GTR (`PartitionCreated.gtr`) in addition to the partition that owns the canonical GTR. When a pipeline mutates the partition's GTR in place after creation but serializes the standalone clone, the written model is stale.

## Active defect: optimize

`optimize --gtr=infer` (the default model) writes an un-normalized `mu` into `gtr.json`.

`packages/treetime/src/optimize/pipeline.rs` captures `let gtr = created.gtr;` before `normalize_partition_rates` mutates the partition's GTR (`gtr.mu /= total_average`, branch lengths scaled by `total_average`). The pipeline returns the pre-normalization clone as `OptimizeOutput.gtr`, while the tree it ships carries the rate-scaled (normalized) branch lengths. Model and tree are mutually inconsistent, and both diverge from v0, which serialized the GTR read back from the partition after normalization (base behavior at `packages/treetime/src/commands/optimize/run.rs` pre-refactor: `write_gtr_json(&partition.read_arc().gtr, ...)` after `normalize_partition_rates`).

For a single partition `total_average == mu`, so after normalization the partition's `gtr.mu == 1.0`. The stale snapshot instead reports the raw inferred `mu`.

## Latent instance: timetree

`packages/treetime/src/timetree/pipeline.rs` (`initialize_partitions_from_params`) also routes `created.gtr` to `TimetreeOutput.gtr`. A scan of `packages/treetime/src/timetree/` finds no `normalize_partition_rates`, `refine_gtr`, or `gtr_mut`: timetree never mutates the GTR after creation, so the snapshot currently equals the live model and output is correct today. It is a latent instance of the same pattern, masked only by timetree not touching the GTR. A future GTR-refinement or normalization step in timetree would silently reopen the defect.

## Correct pattern (already followed by ancestral and prune)

`ancestral` reads `partitions[0].read_arc().gtr().clone()` after `refine_gtr_iterative` (`packages/treetime/src/ancestral/pipeline.rs`); `prune` reads `partitions[0].read_arc().gtr.clone()` (`packages/treetime/src/prune/pipeline.rs`). Both derive the output GTR from the owning partition at output time, never from a standalone clone.

## Root cause

Duplicate source of truth: `PartitionCreated.gtr` is a GTR independent of the partition's owned copy. The two start equal and silently diverge when the partition's copy is mutated. Removing the field and reading the GTR from the owning partition at output time eliminates the divergence and makes the bug unrepresentable.

## Detection gap

The optimize golden master (`packages/treetime/src/optimize/__tests__/test_gm_optimize.rs`) compares only branch lengths and uses a fixed JC69 model; it never exercises the pipeline's output GTR or the `--gtr=infer` normalization, so the regression shipped without a failing test.

## Related issues

- [M-core-partition-init-orchestration-duplication.md](M-core-partition-init-orchestration-duplication.md) -- same `create_marginal_partition` consolidation
- [M-optimize-gtr-mu-coordinate-mismatch.md](M-optimize-gtr-mu-coordinate-mismatch.md) -- separate `mu != 1` scale concern in the optimizer
