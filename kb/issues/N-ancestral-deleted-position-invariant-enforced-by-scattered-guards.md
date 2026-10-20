# Deleted-position invariant is enforced by scattered guards, not by construction

Ancestral reconstruction relies on the invariant that a deleted position carries no character state: every range in `gaps` is also a range in `non_char` (`gaps` is a subset of `non_char`). Fitch leaves (`fn FitchNodeData::new()`) and dense leaves (`fn DenseLeafObs::new()`) hold it at construction [packages/treetime/src/partition/storage/sparse.rs#L248-L266](../../packages/treetime/src/partition/storage/sparse.rs#L248-L266), [packages/treetime/src/partition/storage/dense.rs#L31-L41](../../packages/treetime/src/partition/storage/dense.rs#L31-L41). Internal Fitch nodes did not, so the forward pass could emit a substitution inside a node's own deletion (a position reported as both mutated and gapped on one edge).

The invariant now holds, but each downstream reader re-asserts it, across two representations and both passes, instead of one point that makes a character at a deleted position unrepresentable. The behavior is correct. The concern is that four parallel structures (`sequence`, `non_char`, `gaps`, and the `fitch.variable` / `profile.variable` maps) stay mutually consistent only by pairwise synchronization at several sites, so a future code path that populates `variable` or `profile.variable` at a deleted position without a matching guard silently reintroduces the defect.

## Evidence

Root-cause fix (single point):

- `fn resolve_indels_backward()` runs ahead of sequence construction, and the backward pass unions its resolved gaps into `non_char` where the node's masking ranges are computed [packages/treetime/src/partition/fitch/passes.rs#L130-L137](../../packages/treetime/src/partition/fitch/passes.rs#L130-L137).

Invariant re-asserted at downstream readers (the scatter):

- The forward pass re-unions `non_char` with `gaps` after gap widening [packages/treetime/src/partition/fitch/passes.rs#L205-L214](../../packages/treetime/src/partition/fitch/passes.rs#L205-L214). Correctness depends on running after `fn resolve_indels_forward()`, whose insertion detection reads `non_char` as it stood during the backward pass. Nothing in the code states this order dependency.
- The Fitch substitution resolver skips positions in `non_char` (backward) and in the edge's `gaps` (forward) at three points [packages/treetime/src/partition/fitch/sub.rs#L29](../../packages/treetime/src/partition/fitch/sub.rs#L29), [packages/treetime/src/partition/fitch/sub.rs#L143](../../packages/treetime/src/partition/fitch/sub.rs#L143), [packages/treetime/src/partition/fitch/sub.rs#L165](../../packages/treetime/src/partition/fitch/sub.rs#L165). The backward guard is needed because `variable` is derived from child positions independently of this node's `non_char`.
- Sparse marginal substitution emission in `fn compute_ml_subs_for_nodes()` filters positions that are `non_char` at either endpoint, for parity with the dense `fn edge_subs()` [packages/treetime/src/partition/marginal/sparse/forward.rs#L246](../../packages/treetime/src/partition/marginal/sparse/forward.rs#L246). The stored `profile.variable` still holds an entry whose argmax is an ordinary state at a deleted site; the filter suppresses it only at emit time.
- Sparse marginal reconstruction skips gap positions in `fn variable_states()` before it writes profile states into the sequence [packages/treetime/src/partition/marginal/sparse/reconstruct.rs#L67-L78](../../packages/treetime/src/partition/marginal/sparse/reconstruct.rs#L67-L78), which `fn map_seq_sampled()` and `fn map_seq_overlay()` use [packages/treetime/src/partition/marginal/sparse/reconstruct.rs#L41-L65](../../packages/treetime/src/partition/marginal/sparse/reconstruct.rs#L41-L65). The reconstructed sequence agrees with the indel track, but this is one more read-time filter over a `profile.variable` that contains the deleted site.

## Design concerns

- Not correct by construction. The `fitch.variable` and `profile.variable` structures permit a character state at a deleted position; consistency is a cross-structure property maintained by convention. Any new consumer must know to re-filter.
- Temporal coupling. The forward-pass `non_char` re-union has an undocumented order dependency on `fn resolve_indels_forward()` and its insertion detection.
- Test asserts the symptom, not the invariant. `test_fitch_gap_sub_conflict.rs` checks "no substitution inside a deletion" at the edge output; nothing asserts the subset relation on internal-node structures at construction.

## Mitigating context

The dense path already filters deleted positions at `fn edge_subs()` read time [packages/treetime/src/partition/marginal/dense/partition.rs#L80-L88](../../packages/treetime/src/partition/marginal/dense/partition.rs#L80-L88), so the sparse marginal filter matches an established convention rather than adding a new pattern. The forward Fitch guard on the edge's `gaps` addresses edge-local deletions (parent-inherited gaps that widen on the edge), a distinct condition from node `non_char`, so it is not pure redundancy with the backward union.

## Options

> [!IMPORTANT]
> **Decision required.** The options below differ in scope and in whether the reconstruction representation changes. O1 is correct by construction but changes the profile and indel representation; O2 and O3 keep the current structures. No option is approved.

- O1. Enforce the subset relation at construction of the sparse node and edge data: when a range enters `gaps`, remove it from the `fitch.variable` and `profile.variable` domains in the same operation, so no reader can observe a state at a deleted position. Removes the guards in `fitch/sub.rs`, the filter in `marginal/sparse/forward.rs`, and the gap filter in `fn variable_states()`. Largest change; touches the profile and indel representation.
- O2. Add a single validated combine step that produces `(sequence, non_char, gaps, variable)` as one consistent value with a debug assertion of the subset relation, and have both passes call it. Removes the forward-pass order dependency by making order internal to the combinator. Medium change; keeps the current structures.
- O3. Keep the current guards but add a construction-time invariant check (debug assertion or constructor validation) on internal Fitch nodes and an edge-level assertion that no substitution position lies in the edge's `gaps`. Smallest change; converts the silent-reintroduction risk into a loud failure without redesign.

Behavior must remain unchanged under every option: the reconstruction output and the gap and substitution test results are the contract.

## Validation

- Add a construction-time assertion that every `gaps` range is contained in `non_char` for internal Fitch nodes, and an edge-level assertion that no reported substitution position lies in the edge's deletions. Both must hold on the datasets that motivated the fix.
- Fitch compression, marginal reconstruction (dense and sparse), and the `edge_subs` and `edge_indels` outputs are unchanged.
- No guard removal changes reconstruction output, confirming the guards were enforcing an invariant that construction can own instead.

## Source

Introduced by neherlab/treetime PR #890 "fix: the fitch/parsimony sequence reconstruction did not resolve indels and substitutions correctly", fix commit [`6774b77`](https://github.com/neherlab/treetime/commit/6774b77092e42e31e77628eec102561dba1057aa).

## Related issues

- [N-amino-acid-mutation-indel-representation-undecided.md](N-amino-acid-mutation-indel-representation-undecided.md)
