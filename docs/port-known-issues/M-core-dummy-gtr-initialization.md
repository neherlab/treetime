# Dummy GTR initialization pattern across commands

Three commands (`ancestral`, `optimize`, `prune`) create partition structs with a dummy JC69 GTR model, populate the partition with sequence data, then replace the GTR with the real model. The code is annotated with FIXME comments acknowledging the pattern is undesirable.

## Affected commands

### ancestral

`packages/treetime/src/commands/ancestral/run.rs` creates dummy JC69 for both sparse (line 113) and dense (line 150) partitions, then replaces after `compress_sequences` / `initialize_marginal`:

- [packages/treetime/src/commands/ancestral/run.rs#L113](../../packages/treetime/src/commands/ancestral/run.rs#L113): sparse dummy GTR
- [packages/treetime/src/commands/ancestral/run.rs#L126-L127](../../packages/treetime/src/commands/ancestral/run.rs#L126-L127): sparse GTR replacement
- [packages/treetime/src/commands/ancestral/run.rs#L150](../../packages/treetime/src/commands/ancestral/run.rs#L150): dense dummy GTR
- [packages/treetime/src/commands/ancestral/run.rs#L163-L164](../../packages/treetime/src/commands/ancestral/run.rs#L163-L164): dense GTR replacement

### optimize

`packages/treetime/src/commands/optimize/run.rs` creates dummy JC69 for both sparse (line 63) and dense (line 90) partitions:

- [packages/treetime/src/commands/optimize/run.rs#L63](../../packages/treetime/src/commands/optimize/run.rs#L63): sparse dummy GTR
- [packages/treetime/src/commands/optimize/run.rs#L75-L76](../../packages/treetime/src/commands/optimize/run.rs#L75-L76): sparse GTR replacement
- [packages/treetime/src/commands/optimize/run.rs#L90](../../packages/treetime/src/commands/optimize/run.rs#L90): dense dummy GTR
- [packages/treetime/src/commands/optimize/run.rs#L107-L108](../../packages/treetime/src/commands/optimize/run.rs#L107-L108): dense GTR replacement

### prune

`packages/treetime/src/commands/prune/run.rs` creates dummy JC69 for sparse only (line 64):

- [packages/treetime/src/commands/prune/run.rs#L64](../../packages/treetime/src/commands/prune/run.rs#L64): sparse dummy GTR
- [packages/treetime/src/commands/prune/run.rs#L73-L74](../../packages/treetime/src/commands/prune/run.rs#L73-L74): sparse GTR replacement

## Not affected

The `timetree` command (`commands/timetree/initialization.rs:88-96`) handles this more cleanly: it resolves the GTR model (or creates a JC69 placeholder for `--model=infer`) before constructing partitions, passing the resolved model to the partition constructor. No dummy-then-replace pattern.

## Root cause

`PartitionMarginalSparse` and `PartitionMarginalDense` require a `gtr` field at construction time. GTR inference (`get_gtr_sparse`, `get_gtr_dense`) requires populated partition data (Fitch counts for sparse, marginal profiles for dense). This creates a circular dependency: partitions need a GTR, but GTR inference needs partitions.

The current workaround initializes partitions with a JC69 dummy, populates them, infers the real GTR from the populated data, then replaces the dummy. For `--model=jc69` or other named models, the replacement is a no-op (the dummy could be the real model from the start), but the code path is the same.

## Impact

- For `--model=jc69`: none (dummy equals real model).
- For named models (`hky85`, `jtt92`, etc.): the Fitch parsimony step and `compress_sequences` do not use the GTR model, so the dummy has no effect on sparse partitions. Dense `initialize_marginal` uses the GTR for transition probabilities, so the first marginal pass uses wrong transition matrices. This is corrected by the subsequent `update_marginal` after GTR replacement.
- For `--model=infer`: the inferred model depends on data populated under the dummy. Sparse inference uses Fitch counts (GTR-independent). Dense inference uses profiles from `initialize_marginal` + `update_marginal` under dummy JC69, producing a GTR biased toward equal frequencies. The optimize command mitigates this with an extra `update_marginal` after GTR replacement (added in PR #484 follow-up).

## Proposed solutions

### S1: Deferred GTR field

Make `gtr` an `Option<Gtr>` on partition structs. Construct partitions without a GTR, populate them, infer GTR, then set the field. Operations that need the GTR (`update_marginal`, `expQt`) return an error if GTR is `None`. This makes the dependency explicit in the type system but requires `Option` unwrapping throughout the marginal codepath.

### S2: Two-phase partition construction

Split partition construction into `PartitionBuilder` (no GTR, holds sequence data) and `Partition` (has GTR, ready for inference). The builder phase handles `compress_sequences` and `initialize_marginal` without needing a GTR. The builder converts to a full partition when a GTR is provided. This is a larger refactoring but eliminates the circular dependency structurally.

### S3: Extract GTR inference inputs from raw data

Compute GTR inference inputs (Fitch mutation counts, initial profiles) without constructing a full partition. Pass raw alignment + tree to the GTR inference functions, then construct partitions with the real GTR from the start. This decouples GTR inference from the partition lifecycle but requires extracting the relevant computation from the partition methods.

## v0 handling

v0 has a similar pattern but less visible: `TreeTime.__init__()` accepts a `gtr` parameter (default `'infer'`), and `run()` calls `infer_gtr()` after initial ancestral reconstruction. The GTR is a mutable attribute on the `TreeAnc` object, replaced in place. The circular dependency exists but is hidden by Python's mutable object model.
