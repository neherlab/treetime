# Sparse seq.composition is cached derived state that can go stale

`SparseSeqInfo.composition` is a cached character count derived from `seq.sequence`. It is set during Fitch compression and recomputed during `reconstruct_node_sequence()`, but there is no type-level or runtime enforcement that callers see the correct version. Callers must know whether reconstruction has been called, but nothing in the API prevents reading stale data.

## Problem

`seq.composition` is consumed by:

- GTR inference (`get_mutation_counts_sparse()` at [packages/treetime/src/gtr/infer_gtr/sparse.rs#L51-L71](../../packages/treetime/src/gtr/infer_gtr/sparse.rs#L51-L71)) for Ti (time-in-state) and root_state
- Marginal backward pass leaf seeding at [packages/treetime/src/representation/partition/marginal_passes.rs#L49](../../packages/treetime/src/representation/partition/marginal_passes.rs#L49)
- Marginal backward pass internal `combine_messages()` at [packages/treetime/src/representation/partition/marginal_passes.rs#L88](../../packages/treetime/src/representation/partition/marginal_passes.rs#L88)
- Marginal forward pass `combine_messages()` at [packages/treetime/src/representation/partition/marginal_passes.rs#L194](../../packages/treetime/src/representation/partition/marginal_passes.rs#L194)
- Marginal forward pass child message `fixed_counts` at [packages/treetime/src/representation/partition/marginal_passes.rs#L215](../../packages/treetime/src/representation/partition/marginal_passes.rs#L215)

After `update_marginal()` but before `reconstruct_node_sequence()`:

- `seq.sequence` still holds Fitch-era state (including non-canonical chars like N, R)
- `seq.composition` still holds Fitch-era character counts
- `profile.variable` holds post-marginal MAP posteriors
- `edge_subs_from_graph()` returns MAP-derived mutations (computed dynamically from profiles)

This mixed state means consumers that read `seq.composition` after marginal but before reconstruction get Fitch-era counts while consumers that compute mutations dynamically get MAP-era results. The `get_mutation_counts_sparse()` doc comment documents this prerequisite, but nothing enforces it.

## Current mitigations

- `reconstruct_node_sequence()` recomputes composition from the updated sequence ([packages/treetime/src/representation/partition/marginal_sparse.rs#L473-L477](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L473-L477))
- `get_mutation_counts_sparse()` documents the reconstruction prerequisite
- All production callers currently infer GTR before marginal (Fitch-consistent state), so no mixed provenance occurs in practice
- A `test_reconstruction_changes_composition_for_ambiguous_input` contract test verifies reconstruction causes composition changes for ambiguous sequences

## Why a per-node on-the-fly fix is not possible

`reconstruct_node_sequence()` propagates parent MAP states top-down: a non-root node's sequence starts from the parent's reconstructed sequence, not its own stored sequence. A root ambiguity resolution (N to G) propagates to all descendants at that position even if the descendant has no variable site there. A local per-node composition adjustment cannot replicate this without tree traversal.

## Solutions

### S1. Eliminate stored composition (recommended)

Remove `SparseSeqInfo.composition` entirely. Compute character counts on demand from `seq.sequence` wherever composition is consumed. The marginal passes read composition once per node per pass, and sequences are 1-10KB for typical viral datasets, so the performance cost is negligible.

This eliminates the entire class of staleness bugs: no cached state means no stale state.

**Files affected**: `sparse.rs` (struct definition, constructors), `marginal_passes.rs` (4 reader sites), `sparse.rs` in infer_gtr (2 reader sites), `fitch.rs` (writer sites), `prune/run.rs` (writer site), `marginal_sparse.rs` (reconstruction recomputation).

### S2. Runtime flag

Add `reconstructed: bool` to `SparseNodePartition`. Set `false` during Fitch, set `true` in `reconstruct_node_sequence()`. `get_mutation_counts_sparse()` asserts `reconstructed` when called after marginal (detectable because `profile.variable` is non-empty). Catches misuse at runtime with a clear error message. Does not prevent the staleness, only detects it.

### S3. Type-state pattern

Use different types for pre-reconstruction and post-reconstruction partition states: `PartitionMarginalSparse<Fitch>` vs `PartitionMarginalSparse<Reconstructed>`. GTR inference would only accept `<Reconstructed>`. Catches misuse at compile time but adds generic complexity throughout the partition trait hierarchy.

## v0 comparison

v0 does not have this problem. `node.cseq` (compressed sequence) is overwritten after each marginal pass, and the `mutations` property dynamically compares current `node.up.cseq` vs `node.cseq`. Composition data always reflects the current reconstruction state.

## Related

- [M-gtr-per-site-rate-variation](M-gtr-per-site-rate-variation.md) - per-site rate variation, another GTR inference issue
- [M-core-dummy-gtr-initialization](M-core-dummy-gtr-initialization.md) - dummy GTR pattern that shapes when GTR inference runs relative to marginal
