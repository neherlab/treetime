# Deleted-position invariant is enforced by scattered guards, not by construction

Ancestral reconstruction relies on the invariant that a deleted position carries no character state: every range in `gaps` is also a range in `non_char` (`gaps` is a subset of `non_char`). Leaves (`SparseNodePartition::new`) and the dense representation (`DenseSeqInfo::new`) hold it at construction. Internal sparse nodes did not, so the forward pass could emit a substitution inside a node's own deletion (a position reported as both mutated and gapped on one edge).

The fix restores the invariant, but enforces it by re-asserting it at each downstream reader across two representations and both passes, rather than making a character-at-a-deleted-position unrepresentable at a single point. The behavior is now correct; the concern is that four parallel structures (`sequence`, `non_char`, `gaps`, and the `variable` / `profile.variable` maps) are kept mutually consistent by pairwise synchronization at several sites, so a future code path that populates `variable`/`profile` at a deleted position without a matching guard silently reintroduces the defect.

## Evidence

Root-cause fix (correct, single chokepoint):

- Backward pass unions resolved gaps into `non_char` where the node's masking ranges are computed [`packages/treetime/src/ancestral/fitch.rs#L188`](../../packages/treetime/src/ancestral/fitch.rs#L188), after reordering `resolve_indels_backward` ahead of sequence construction [`packages/treetime/src/ancestral/fitch.rs#L178`](../../packages/treetime/src/ancestral/fitch.rs#L178).

Invariant re-asserted at downstream readers (the scatter):

- Forward pass re-unions `non_char` with `gaps` after gap-widening; correctness depends on running after `resolve_indels_forward`, whose insertion detection reads `non_char` as it stood during the backward pass [`packages/treetime/src/ancestral/fitch.rs#L288`](../../packages/treetime/src/ancestral/fitch.rs#L288). This is call-order coupling justified by comment.
- Fitch substitution resolver skips positions in `non_char` (backward) and in the edge's `gaps` (forward) at three points [`packages/treetime/src/ancestral/fitch_sub.rs#L33`](../../packages/treetime/src/ancestral/fitch_sub.rs#L33), [`packages/treetime/src/ancestral/fitch_sub.rs#L150`](../../packages/treetime/src/ancestral/fitch_sub.rs#L150), [`packages/treetime/src/ancestral/fitch_sub.rs#L173`](../../packages/treetime/src/ancestral/fitch_sub.rs#L173). The backward guard is needed because `variable` is derived from child positions independently of this node's `non_char`.
- Sparse marginal substitution emission filters positions that are `non_char` at either endpoint, for parity with the dense `edge_subs()` [`packages/treetime/src/partition/marginal_passes.rs#L246`](../../packages/treetime/src/partition/marginal_passes.rs#L246). The stored `profile.variable` still holds a residue whose argmax is an ordinary state at a deleted site; the filter suppresses it only at emit time.
- Sparse marginal reconstruction writes residues from the profile and then re-fills deletion ranges with gaps to overwrite them [`packages/treetime/src/partition/marginal_sparse.rs#L116`](../../packages/treetime/src/partition/marginal_sparse.rs#L116). This repairs the output after producing the inconsistent state rather than not producing it.

## Design concerns

- Not correct by construction. The `variable` and `profile.variable` structures permit a character state at a deleted position; consistency is a cross-structure property maintained by convention. Any new consumer must know to re-filter.
- Temporal coupling. The forward-pass `non_char` re-union has an order dependency on `resolve_indels_forward` and insertion detection, documented in a comment instead of removed.
- Output repair over prevention. `reconstruct_map_seq_sampled` writes then overwrites deleted sites, so the profile-derived sequence transiently disagrees with the indel track.
- Test asserts the symptom, not the invariant. `test_fitch_gap_sub_conflict.rs` checks "no substitution inside a deletion" at the edge output; nothing asserts the subset relation on internal-node structures at construction.

## Mitigating context

The dense path already filters deleted positions at `edge_subs()` read time, so the sparse marginal filter [`packages/treetime/src/partition/marginal_passes.rs#L246`](../../packages/treetime/src/partition/marginal_passes.rs#L246) matches an established convention rather than adding a new pattern. The forward `fitch_sub` guard on the edge's `gaps` addresses edge-local deletions (parent-inherited gaps that widen on the edge), a genuinely distinct condition from node `non_char`, so it is not pure redundancy with the backward union.

## Options

- O1. Enforce the subset relation at construction of the sparse node/edge partition: when a range enters `gaps`, remove it from the `variable` and `profile.variable` domains in the same operation, so no reader can observe a residue at a deleted position. Removes the `fitch_sub` and `marginal_passes` guards and the `marginal_sparse` overwrite. Largest change; touches the profile/indel representation.
- O2. Add a single validated combine step that produces `(sequence, non_char, gaps, variable)` as one consistent value with a debug assertion of the subset relation, and have both passes call it. Removes the forward-pass ordering comment by making order internal to the combinator. Medium change; keeps the current structures.
- O3. Keep the current guards but add a construction-time invariant check (debug assertion or constructor validation) on internal sparse nodes and an edge-level assertion that no substitution position lies in the edge's `gaps`. Smallest change; converts the silent-reintroduction risk into a loud failure without redesign.
- O4. Replace the `marginal_sparse` write-then-overwrite with not writing residues at deleted positions in the profile-application loop. Local improvement, independent of O1-O3.

O3 and O4 are compatible with the current design and reduce the latent risk immediately. O1 is the correct-by-construction target and needs approval because it changes the reconstruction representation. Behavior must remain unchanged: the reconstruction output and the gap/substitution test results are the contract.

## Validation

- Add a construction-time assertion that every `gaps` range is contained in `non_char` for internal sparse nodes, and an edge-level assertion that no reported substitution position lies in the edge's deletions. Both must hold on the datasets that motivated the fix.
- Fitch compression, marginal reconstruction (dense and sparse), and `edge_subs`/`edge_indels` outputs are unchanged versus the current branch.
- No guard removal changes reconstruction output, confirming the guards were enforcing an invariant that construction can own instead.

## Source

Introduced by neherlab/treetime PR #890 "fix: the fitch/parsimony sequence reconstruction did not resolve indels and substitutions correctly", fix commit [`6774b77`](https://github.com/neherlab/treetime/commit/6774b77092e42e31e77628eec102561dba1057aa).

## Related issues

- [M-partition-compressed-exposes-fitch-storage.md](M-partition-compressed-exposes-fitch-storage.md)
- [H-gtr-model-state-can-break-invariants.md](H-gtr-model-state-can-break-invariants.md)
- [N-amino-acid-mutation-indel-representation-undecided.md](N-amino-acid-mutation-indel-representation-undecided.md)
