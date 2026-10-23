# Sparse edge merges compose substitutions and indels separately

When a topology operation merges two sparse edges into one, it composes the substitutions of the two edges and their indels in two independent steps. A site with a substitution on one edge and an indel on the other keeps both events on the merged edge, and the indel keeps the character of the wrong end. The merged edge rebuilds the child correctly in the forward direction, but its inversion rebuilds the wrong state, so a later reroot through the merged edge changes the sequences of nodes it only moves.

## Reproduction

Tree `(((Y1:0.1,Y2:0.1)Y:0.1,Z:0.1)X:0.1,W:0.1,V:0.1)root;`, sequences `Y1=-CG`, `Y2=-CG`, `Z=GCG`, `W=ACG`, `V=ACT`, sparse Fitch reconstruction through production entry points (`create_fitch_partition()`, `collapse_edge()`, `apply_root_at_edge()`, `PartitionSet::apply_reroot()`):

1. Fitch puts the substitution `A1G` on root->X and the deletion `0--1: G -> -` on X->Y
2. `collapse_edge()` on root->X gives root->Y the substitution `A1G` and the deletion `0--1: G -> -`. In the forward direction the substitution sets `G` and the deletion removes it, so Y rebuilds as `-CG` (correct)
3. A reroot at Y->Y1 (split 0.5) inverts root->Y to the substitution `G1A` and the insertion of `G`. The rebuild of the old root applies the substitution (`A`), then the insertion (`G`)

| Node | Before the reroot | After the reroot |
| ---- | ----------------- | ---------------- |
| root | `ACG` | `GCG` |
| W | `ACG` | `GCG` |
| V | `ACT` | `GCT` |

A merged edge that keeps both events records the deletion with the state `G` of the intermediate node X, where the parent root has `A`; the composed edge from root to Y is a deletion of `A` without a substitution.

The same separation affects an insertion followed by a substitution at its site (parent edge: insertion of `A`; child edge: substitution `A->G`): the merged edge stores both, and `parsimony_seq()` applies the substitution first and the insertion second, so the child rebuilds as `A` instead of `G`.

## Affected code

- `fn collapse_edge()` [packages/treetime/src/optimize/topology/collapse.rs#L28-L37](../../packages/treetime/src/optimize/topology/collapse.rs#L28-L37) chains the substitutions with `chain_fitch_subs()` and composes the indels with `composed_edge_pieces()` separately
- `fn merge_edge_pair()` [packages/treetime/src/partition/marginal/sparse/alignment_indels.rs#L76-L96](../../packages/treetime/src/partition/marginal/sparse/alignment_indels.rs#L76-L96) does the same for the trivial old root of a reroot and for `PartitionSet::reconcile_topology()`
- `fn compose_substitutions()` [packages/treetime/src/seq/mutation.rs#L305](../../packages/treetime/src/seq/mutation.rs#L305) and `fn compose_indels()` [packages/treetime/src/seq/indel.rs#L9](../../packages/treetime/src/seq/indel.rs#L9) each see only one event kind
- `fn parsimony_seq()` [packages/treetime/src/partition/marginal/sparse/reconstruct.rs#L12-L37](../../packages/treetime/src/partition/marginal/sparse/reconstruct.rs#L12-L37) applies all substitutions, then all indels; `fn apply_edge_to_sequence()` in `packages/treetime/src/partition/marginal/sparse/reroot.rs` undoes them in the same order

## Frequency

`collapse_edge()` runs for `prune --collapse-short`, `--collapse-empty` and named nodes [packages/treetime/src/prune/prune.rs#L76](../../packages/treetime/src/prune/prune.rs#L76), for edges whose optimal length is zero in the optimize loop [packages/treetime/src/optimize/run_loop.rs#L380](../../packages/treetime/src/optimize/run_loop.rs#L380), and for polytomy resolution [packages/treetime/src/optimize/topology/resolve_polytomy.rs#L195](../../packages/treetime/src/optimize/topology/resolve_polytomy.rs#L195). The reproduction needs a short edge with a substitution above an edge with a deletion at the same site, followed by a reroot through the merged edge.

> [!IMPORTANT]
> **Investigation required.** Collect whether Fitch and the topology operations produce the insertion-then-substitution pair in `merge_edge_pair()` at a reroot: a bifurcating Fitch root takes the state of its child without a gap, so a substitution on root->Y at a site where X has a gap needs an earlier operation. Measure how often the sub-then-deletion pair occurs on the example datasets with `prune --collapse-short` or `optimize` followed by a reroot.

## Fix direction

Compose the events of a merged edge site by site over all event kinds (substitutions, indels and the unknown ranges of the removed node), so that the merged edge holds the net change from its parent to its child. The fix belongs to one composition that every sparse topology operation (collapse, trivial root merge, reroot inversion, reversion hoist) uses; see also [N-sparse-inverted-indel-pair-rebuilds-wrong-state.md](N-sparse-inverted-indel-pair-rebuilds-wrong-state.md).
