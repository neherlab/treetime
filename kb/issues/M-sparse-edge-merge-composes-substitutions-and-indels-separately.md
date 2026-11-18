# Sparse edge merges compose substitutions and indels separately

When a topology operation merges two sparse edges into one, it composes the substitutions of the two edges and their indels in two independent steps. A site with a substitution on one edge and an indel on the other keeps the indel with the character of the wrong end, and keeps both events where the substitution survives the merge. The merged edge rebuilds the child correctly in the forward direction, but its inversion rebuilds the wrong state, so a later reroot through the merged edge changes the sequences of nodes it only moves.

## Reproduction

Tree `(((Y1:0.1,Y2:0.1)Y:0.1,Z:0.1)X:0.1,W:0.1,V:0.1)root;`, sequences `Y1=-CG`, `Y2=-CG`, `Z=GCG`, `W=ACG`, `V=ACT`, sparse Fitch reconstruction through production entry points (`create_fitch_partition()`, `collapse_edge()`, `apply_root_at_edge()`, `PartitionSet::apply_reroot()`):

1. Fitch puts the substitution `A1G` on root->X and the deletion `0--1: G -> -` on X->Y
2. `collapse_edge()` on root->X drops the composed substitution `A1G`, because Y has no data at site 1, and gives root->Y the deletion `0--1: G -> -`. In the forward direction the deletion removes the character, so Y rebuilds as `-CG` (correct)
3. A reroot at Y->Y1 (split 0.5) inverts root->Y to the insertion of `G`, so the rebuild of the old root puts `G` at site 1. `merge_edge_pair()` keeps both events of such a pair; the inversion then gives the substitution `G1A` and the insertion of `G`, and the rebuild applies the substitution (`A`), then the insertion (`G`), with the same result

| Node | Before the reroot | After the reroot |
| ---- | ----------------- | ---------------- |
| root | `ACG` | `GCG` |
| W | `ACG` | `GCG` |
| V | `ACT` | `GCT` |

The merged edge records the deletion with the state `G` of the intermediate node X, where the parent root has `A`; the composed edge from root to Y is a deletion of `A` without a substitution.

The same separation affects an insertion followed by a substitution at its site (parent edge: insertion of `A`; child edge: substitution `A->G`): the merged edge stores both, and `parsimony_seq()` applies the substitution first and the insertion second, so the child rebuilds as `A` instead of `G`.

## Affected code

- `fn collapse_edge()` [packages/treetime/src/optimize/topology/collapse.rs#L43-L70](../../packages/treetime/src/optimize/topology/collapse.rs#L43-L70) chains the substitutions with `chain_fitch_subs()`, drops those at sites where the child has no data, and composes the indels with `composed_edge_pieces()` separately
- `fn merge_edge_pair()` [packages/treetime/src/partition/marginal/sparse/alignment_indels.rs#L48-L68](../../packages/treetime/src/partition/marginal/sparse/alignment_indels.rs#L48-L68) does the same for the trivial old root of a reroot and for `PartitionSet::reconcile_topology()`
- `fn compose_substitutions()` [packages/treetime/src/seq/mutation.rs#L345](../../packages/treetime/src/seq/mutation.rs#L345) and `fn compose_indels()` [packages/treetime/src/seq/indel.rs#L9](../../packages/treetime/src/seq/indel.rs#L9) each see only one event kind
- `fn parsimony_seq()` [packages/treetime/src/partition/marginal/sparse/reconstruct.rs#L23-L48](../../packages/treetime/src/partition/marginal/sparse/reconstruct.rs#L23-L48) applies all substitutions, then all indels; `fn apply_edge_to_sequence()` in `packages/treetime/src/partition/marginal/sparse/reroot.rs` undoes them in the same order

## Frequency

`collapse_edge()` runs for `prune --collapse-short`, `--collapse-empty` and named nodes [packages/treetime/src/prune/prune.rs#L77](../../packages/treetime/src/prune/prune.rs#L77), for edges whose optimal length is zero in the optimize loop [packages/treetime/src/optimize/run_loop.rs#L329](../../packages/treetime/src/optimize/run_loop.rs#L329), and for polytomy resolution and reversion hoists [packages/treetime/src/optimize/topology/resolve_polytomy.rs#L136](../../packages/treetime/src/optimize/topology/resolve_polytomy.rs#L136), [packages/treetime/src/optimize/topology/resolve_polytomy.rs#L225](../../packages/treetime/src/optimize/topology/resolve_polytomy.rs#L225). The reproduction needs a short edge with a substitution above an edge with a deletion at the same site, followed by a reroot through the merged edge.

> [!IMPORTANT]
> **Investigation required.** Collect whether Fitch and the topology operations produce the insertion-then-substitution pair in `merge_edge_pair()` at a reroot: a bifurcating Fitch root takes the state of its child without a gap, so a substitution on root->Y at a site where X has a gap needs an earlier operation. Measure how often the sub-then-deletion pair occurs on the example datasets with `prune --collapse-short` or `optimize` followed by a reroot.

## Fix direction

Compose the events of a merged edge site by site over all event kinds (substitutions, indels and the unknown ranges of the removed node), so that the merged edge holds the net change from its parent to its child. The fix belongs to one composition that every sparse topology operation (collapse, trivial root merge, reroot inversion, reversion hoist) uses; see also [N-sparse-inverted-indel-pair-rebuilds-wrong-state.md](N-sparse-inverted-indel-pair-rebuilds-wrong-state.md).

## Proposed implementation plan

1. Add one site-wise composition of two sparse edges, for example `fn compose_edge_events(parent: &SparseEdgeObs, child: &SparseEdgeObs, parent_seq_at, ...) -> SparseEdgeObs` beside `compose_indels()` in `packages/treetime/src/seq/`: for each site touched by either edge, compute the state of the parent of the merged edge and the state of its child, and emit one net event (none, a substitution, a deletion with the parent's character, or an insertion with the child's character), together with the unknown ranges of the removed node
2. Use it in `fn collapse_edge()` (`packages/treetime/src/optimize/topology/collapse.rs`) and in `fn merge_edge_pair()` (`packages/treetime/src/partition/marginal/sparse/alignment_indels.rs`), in place of the pair `chain_fitch_subs()` + `composed_edge_pieces()` and of the no-data filter of `collapse_edge()`, so that every sparse topology operation (collapse, trivial root merge, reconcile, reversion hoist) shares it
3. Keep `parsimony_seq()` and `apply_edge_to_sequence()` as they are; after step 1 a site carries at most one event, so their order no longer matters. Optionally assert this invariant in debug builds
4. The marginal estimates (`compose_substitutions()` on `reconstruction.edges.estimates`) need the same rule when an estimate site meets an indel of the other edge; check this in the same change

Tests (in `packages/treetime/src/optimize/topology/__tests__/` and `packages/treetime/src/partition/marginal/sparse/__tests__/`):

- Regression test of the reproduction above: after `collapse_edge()` and `apply_root_at_edge()`, the sequences of root, W and V stay `ACG`, `ACG`, `ACT`. Oracle: the input sequences (a reroot moves nodes, it does not change them)
- Unit cases of the composition (rstest table): substitution then deletion, insertion then substitution, deletion then insertion, substitution then substitution back to the parent state. Oracle: hand-derived net change from parent to child
- Property test: for random sparse trees, merging an edge and then rebuilding with `parsimony_seq()` gives the same child sequence as the two edges applied in order, and the inverted merged edge rebuilds the parent. Oracle: the unmerged path

Validation: `just test-rs`, `just lint-rs`, then `just smoke` on `prune --collapse-short` and `optimize` per `cli-verification.instructions.md`. The open investigation (how often the pairs occur) does not block the fix.
