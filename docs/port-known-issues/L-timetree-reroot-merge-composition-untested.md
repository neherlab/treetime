# Reroot merge branch in `PartitionMarginalSparse::apply_reroot` lacks direct assertion

The merge branch of `PartitionMarginalSparse::apply_reroot` at [packages/treetime/src/representation/partition/marginal_sparse.rs#L480-L506](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L480-L506) runs when rerooting removes a trivial old root. It composes substitutions via `compose_substitutions` and concatenates indel lists into the merged edge.

## Current coverage

- `test_sparse_reroot_split_root_sequence_applies_path_node_masks` at [packages/treetime/src/commands/timetree/optimization/**tests**/test_reroot.rs#L163](../../packages/treetime/src/commands/timetree/optimization/__tests__/test_reroot.rs#L163): split-root sequence masking
- `test_sparse_reroot_inverts_subs_and_indels_on_path`: sub inversion and indel inversion on the reroot path
- `test_reroot_tree_sparse_with_edge_split`: end-to-end reroot smoke path

## What is not directly asserted

None of the existing tests:

- Force a trivial old root that requires merge branch execution
- Assert the composed substitution list on the merged edge matches the concatenation-then-inversion of the two constituent edges' subs
- Assert the indel list on the merged edge matches the concatenation of the constituent edges' indels
- Assert the merged edge carries no stale metadata from either constituent edge

## Latent correctness risk: overlapping indels

The merged edge stores a concatenated indel list. `edge_state_from_parent` at [packages/treetime/src/representation/partition/marginal_sparse.rs#L161-L177](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L161-L177) selects a covering indel with `Iterator::find`, returning the first match rather than composing all matches in order. When both constituent edges carry indels spanning the same aligned position (for example, an insertion on the parent-side edge that the child-side edge then deletes), only the parent-side indel is observed. The merged edge's apparent child state at the overlapping positions is the parent-side transition, not the composed parent-then-child transition.

### Decision required

The fix is not uniquely determined by the observation. An implementer must choose one of:

- **O1. Enforce non-overlap at merge time**: extend the merge branch of `apply_reroot` to canonicalize the concatenated indel list into non-overlapping intervals (insert-then-delete collapses; adjacent same-type indels merge). `edge_state_from_parent` stays with `.find`.
- **O2. Compose on lookup**: keep the merge branch as concatenation; change `edge_state_from_parent` to fold all covering indels in list order so the last one wins at overlapping positions.
- **O3. Both**: canonicalize on merge AND make `edge_state_from_parent` robust to overlap for defense in depth.

The precedence contract the lookup must match lives in `child_canonical_state_from_parent_state` at [packages/treetime/src/representation/partition/marginal_sparse.rs#L215](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L215) with a parameterized test at [packages/treetime/src/commands/ancestral/**tests**/test_marginal_sparse.rs#L642](../../packages/treetime/src/commands/ancestral/__tests__/test_marginal_sparse.rs#L642). That helper handles single-edge precedence (gap > unknown > sub > parent_state); the merge-branch decision is about multi-indel composition on a single edge, which the helper does not currently cover.

## How to trigger the merge branch from a test

The merge runs when `remove_trivial_root` is set on `RerootParams`. `RerootParams` is defined at [packages/treetime/src/commands/clock/reroot.rs#L17](../../packages/treetime/src/commands/clock/reroot.rs#L17) with `remove_trivial_root` at [`reroot.rs#L26`](../../packages/treetime/src/commands/clock/reroot.rs#L26). The condition that emits `edge_merge` into `RerootResult` is at [`reroot.rs#L161`](../../packages/treetime/src/commands/clock/reroot.rs#L161).

A test driving the real path:

- Build a tree where the old root has exactly two children and rerooting to one of them leaves the old root with a single remaining child, making it trivial.
- Seed partition data so each relevant edge carries known `subs` and/or `indels` to test composition.
- Call `reroot_in_place` with `RerootParams { remove_trivial_root: true, ... }`, feed the returned `RerootResult.changes` into `apply_reroot`.
- Assert the merged edge's `subs` via `compose_substitutions` semantics and `indels` via the chosen canonical form (O1/O2/O3).

Alternative (bypasses the producer): hand-construct `RerootChanges` with `edge_merge: Some(EdgeMergeInfo { ... })` and call `apply_reroot` directly. This pattern already exists in `test_reroot.rs` for the inversion path and is easier to assemble; the tradeoff is that it does not exercise `reroot_in_place` metadata production.

## Proposed test

A focused reroot integration test that:

1. Constructs a three-node path (root, intermediate, new root candidate) with known substitutions and indels on both edges, including an overlapping indel pair
2. Triggers reroot with `remove_trivial_root = true` so the intermediate node is removed and the two edges merge
3. Asserts the merged edge's `subs` and `indels` match the expected composition against an independent oracle computed by hand
4. Asserts `edge_state_from_parent` through the merged edge reproduces the expected child state, including at overlapping indel positions

## What it catches

- Errors in substitution composition (`compose_substitutions`) during merge
- Errors in indel concatenation across merged edges, including overlapping indels that the current lookup cannot resolve
- Stale edge metadata (messages, transmission) not reset after merge
