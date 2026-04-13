# Sparse reroot integration contract and end-to-end path are not directly asserted

## Cached `RerootResult` metadata is not verified against a graph-path oracle

`reroot_in_place` at [packages/treetime/src/commands/clock/reroot.rs#L112](../../packages/treetime/src/commands/clock/reroot.rs#L112) populates `RerootResult.old_root_key`, `inverted_edge_keys`, `path_node_keys`, `edge_split`, and `edge_merge`. Sparse partition updates consume these fields directly rather than recomputing the reroot path locally.

Existing reroot tests in `packages/treetime/src/commands/clock/__tests__/test_reroot.rs` assert node counts and old-root preservation but do not inspect the cached metadata. Sparse partition tests in `packages/treetime/src/commands/timetree/optimization/__tests__/test_reroot.rs` hand-construct `RerootChanges` values, bypassing the producer.

No test asserts that the cached metadata emitted by `reroot_in_place` matches a fresh graph-path oracle computed from the post-reroot topology.

### Graph-path oracle construction

Independent oracle for the cached metadata:

- **`old_root_key`**: capture before rerooting, assert unchanged by comparing to `graph.get_exactly_one_root()?.read_arc().key()` recorded on the pre-reroot graph.
- **`inverted_edge_keys`**: after rerooting, walk the post-reroot graph from the new root down to the old-root key. At each step, record the edge key. The reverse of that list is the expected set of edges that had to be inverted. Alternatively walk from old root toward new root on a clone of the pre-reroot graph.
- **`path_node_keys`**: the same walk recorded at node granularity, starting at `old_root_key` and ending at `new_root_key`.
- **`edge_split`**: if `RerootParams.split_edge` is set, expect `Some(EdgeSplitInfo)` where `old_edge_key` identifies the edge chosen by `find_best_root`, and `new_node_key` is the freshly inserted node. The test sets up the graph so only one best edge qualifies and asserts on its fields.
- **`edge_merge`**: if `RerootParams.remove_trivial_root` is set and the old root becomes trivial after topology change, expect `Some(EdgeMergeInfo)` where `parent_edge_key` and `child_edge_key` are the two edges incident to the old root and `merged_edge_key` is the new edge replacing them.

The oracle is ordinary graph traversal plus pre/post snapshots. No shared code with `reroot_in_place` itself.

## `reroot_tree` producer-to-consumer handoff is only smoke-tested

`reroot_tree` at [packages/treetime/src/commands/timetree/optimization/reroot.rs#L19](../../packages/treetime/src/commands/timetree/optimization/reroot.rs#L19) forwards cached reroot metadata from clock reroot directly into partition updates. Existing tests either use hand-constructed `RerootChanges` or smoke-test that rerooting does not panic and leaves a plausible tree.

No single test exercises the real handoff where `reroot_in_place` emits metadata, `reroot_tree` consumes it, and post-reroot reconstructed sequences are asserted against a path-walk oracle.

### End-to-end test recipe

1. Build a tree with known sequences at every leaf and a known starting root.
2. Run Fitch compression and one marginal update to populate sparse partitions.
3. Capture reconstructed sequences via `ancestral_reconstruction_marginal`.
4. Call `reroot_in_place` to reroot to a chosen internal or leaf node.
5. Call `reroot_tree` feeding the returned `RerootResult`.
6. Capture reconstructed sequences again after rerooting.
7. Assert the sequences at every node are identical to the pre-reroot sequences up to the orientation change (the canonical state at each position depends only on the root sequence and edge deltas, which reroot inverts consistently).

Assertion oracle: the pre-reroot sequence at each node is the ground truth because canonical state is topology-invariant when edges are inverted correctly.

## `reconstruct_split_root_sequence` is verified only on mask-only paths

`test_sparse_reroot_split_root_sequence_applies_path_node_masks` at [packages/treetime/src/commands/timetree/optimization/**tests**/test_reroot.rs#L163](../../packages/treetime/src/commands/timetree/optimization/__tests__/test_reroot.rs#L163) covers path-node gap and unknown masks but uses a path with no substitutions or insertions on the inverted edges.

`reconstruct_split_root_sequence` at [packages/treetime/src/representation/partition/marginal_sparse.rs#L382](../../packages/treetime/src/representation/partition/marginal_sparse.rs#L382) applies `edge_state_from_parent` per inverted edge and then overlays path-node masks. The substitution-plus-indel path through this loop is not directly asserted.

The per-position precedence the split-root reconstruction must match is the contract tested for the single-edge helper at [packages/treetime/src/commands/ancestral/**tests**/test_marginal_sparse.rs#L642](../../packages/treetime/src/commands/ancestral/__tests__/test_marginal_sparse.rs#L642) (gap > unknown > sub > parent_state). The split-root test should extend the same precedence to multi-edge paths.

### Hand-walk oracle for multi-edge paths

For a path `root -> n1 -> n2 -> ... -> new_root` with edges `e0, e1, ..., ek`:

- Start with `sequence = root.seq.sequence.clone()`.
- For each edge `e_i` in order: for each position `pos`, apply `edge_state_from_parent(e_i, pos, sequence[pos])`. This uses the edge's subs and indels to derive the state at the node on the other side of the edge.
- After each edge is applied, overlay the masks of the intermediate node that the edge just brought us to: fill gap ranges with `alphabet.gap()` and unknown ranges with `alphabet.unknown()`.
- The final `sequence` is the expected reconstructed sequence at `new_root`.

This is exactly what `reconstruct_split_root_sequence` claims to do. The test constructs the expected sequence by the same procedure and asserts equality.

## Timetree sparse contribution wrappers are not directly tested

`collect_edge_contributions` at [packages/treetime/src/commands/timetree/inference/branch_length_likelihood.rs#L14](../../packages/treetime/src/commands/timetree/inference/branch_length_likelihood.rs#L14) and `collect_contributions` at [packages/treetime/src/commands/timetree/inference/runner.rs#L137](../../packages/treetime/src/commands/timetree/inference/runner.rs#L137) own graph-aware sparse contribution assembly for timetree inference. No focused test calls either wrapper.

Underlying `create_edge_contribution` is covered by [packages/treetime/src/commands/optimize/**tests**/test_dense_sparse_equivalence/test_dense_sparse_equivalence_initial.rs#L70](../../packages/treetime/src/commands/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_initial.rs#L70). The wrapper-level invariants (contribution count and ordering, mixed dense+sparse fan-in, per-partition cache isolation) are not asserted.

### Multi-partition fixture pattern

Neither the timetree inference code nor existing tests construct a mixed dense+sparse partition list, but the type `PartitionOptimizeVec = Vec<Arc<RwLock<dyn PartitionOptimizeOps>>>` admits it. A wrapper test can:

- Build `PartitionMarginalDense` from a small alignment with `alphabet = Alphabet::Nuc`.
- Build `PartitionMarginalSparse` from the same alignment.
- Cast each to `Arc<RwLock<dyn PartitionOptimizeOps>>` and push into a `PartitionOptimizeVec` in a known order.
- Call `collect_edge_contributions` or `collect_contributions`.
- Assert the returned contributions include one entry per `(edge, partition)` pair in the expected order, matching direct calls to `create_edge_contribution` for each.

Cache isolation check: capture the `ExactStateCache` state after processing a sparse partition, process a second sparse partition, assert the caches do not share entries. The current test at `test_optimize_method.rs` indirectly exercises isolation but does not assert non-sharing.

## Proposed tests

- Direct `reroot_in_place` test that compares `RerootResult.inverted_edge_keys`, `RerootResult.path_node_keys`, `old_root_key`, `edge_split`, and `edge_merge` against a fresh graph-path oracle computed as described above.
- End-to-end `reroot_tree` test following the recipe above: reconstruct before and after, assert sequences preserved.
- Extended split-root test that builds an inverted-edge path carrying substitutions and indels plus intermediate-node masks, with an expected sequence computed by the hand-walk oracle.
- Wrapper tests that build one dense and one sparse partition, call `collect_edge_contributions` and `collect_contributions`, and compare returned contributions against direct `create_edge_contribution` calls.

## What it catches

- Off-by-one path capture in `reroot_in_place` metadata
- Stale or out-of-order edge references in the cached path
- Wrong split-root sequence from mutation-bearing paths
- Dropped sparse partitions, duplicate contributions, or cross-partition cache leakage in timetree wrappers
