# Inverting a composed deletion and insertion over one range rebuilds the wrong state

`compose_indels()` keeps a deletion of the parent edge and an insertion of the child edge over the same range when their characters differ. The pair rebuilds the child correctly in list order, but the inversion of an edge inverts the kind of each indel and keeps the list order, so the inverted pair rebuilds a gap where the parent has a character.

## Reproduction

Function-level, with the production functions (`compose_indels()`, `InDel::invert()`, `parsimony_seq()`), run once as a temporary test:

1. `compose_indels([0--1: A -> -], [0--1: - -> G])` returns both events: `["0--1: A -> -", "0--1: - -> G"]`
2. Applied to the parent `AC` in list order, the pair rebuilds the child `GC` (correct)
3. Inverted in place (insertion of `A`, then deletion of `G`) and applied to the child `GC`, it rebuilds `-C` where the parent is `AC`

## Affected code

- `fn compose_indels()` [packages/treetime/src/seq/indel.rs#L109-L118](../../packages/treetime/src/seq/indel.rs#L109-L118) keeps both events of a deletion and an insertion over one range unless range and characters are equal
- `fn invert_edge()` [packages/treetime/src/partition/marginal/sparse/reroot.rs#L84-L104](../../packages/treetime/src/partition/marginal/sparse/reroot.rs#L84-L104) and `fn remove_stem_root()` (same file) invert each indel in place and keep the order
- `fn parsimony_seq()` [packages/treetime/src/partition/marginal/sparse/reconstruct.rs#L12-L37](../../packages/treetime/src/partition/marginal/sparse/reconstruct.rs#L12-L37) and `fn apply_edge_to_sequence()` (`packages/treetime/src/partition/marginal/sparse/reroot.rs`) apply the indels in list order

> [!IMPORTANT]
> **Investigation required.** No production path is known to produce the pair. Fitch never stores an insertion, because an internal node with a child that has a character takes a character, so an insertion appears only on an edge that a reroot inverted. The pair then needs a deletion on a parent edge and an insertion on its child edge through the same node, composed by a collapse or a merge after a reroot (for example in the reroot iterations of `timetree`). Collect a reproduction through `apply_root_at_edge()` followed by `collapse_edge()` or `merge_edge_pair()` before the fix, or show that no sequence of operations reaches it.

## Fix direction

The fix belongs to the composition of all event kinds that [M-sparse-edge-merge-composes-substitutions-and-indels-separately.md](M-sparse-edge-merge-composes-substitutions-and-indels-separately.md) describes: a composed edge that holds the net change over each range (here a substitution `A->G` instead of a deletion and an insertion) inverts without depending on list order.
