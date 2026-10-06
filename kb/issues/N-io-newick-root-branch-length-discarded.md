# Root branch length silently discarded during Newick parsing

Severity: N (negligible)
Category: I/O
Crate: `util-newick`

## Description

The Newick parser accepts root-level branch lengths (`(A:0.1,B:0.2):0.5;`) but silently discards the root's `:0.5`. `NewickGraph` has no field to represent root branch lengths since the root node has no incoming edge. Comments after the root's `:` become attributes of the root node.

Round-trip of trees with root branch lengths loses data: `(A:0.1,B:0.2):0.5;` parses and writes back as `(A:0.1,B:0.2);`.

## Impact

Root branch lengths appear in some tool outputs (IQ-TREE, BEAST) and can carry evolutionary distance information (e.g., stem branch above displayed root). The loss is silent with no warning to callers.

## Fix

1. Add `root_branch_length: Option<f64>` to `NewickGraph`
2. Store the root's branch length in `fn Builder::finish()` in [packages/util-newick/src/parse.rs](../../packages/util-newick/src/parse.rs) instead of discarding it
3. Write the stored value after the root subtree in `fn write_newick()` in [packages/util-newick/src/write.rs](../../packages/util-newick/src/write.rs)
4. Carry the value through `fn graph_from_newick()` in [packages/treetime-io/src/nwk.rs](../../packages/treetime-io/src/nwk.rs) and decide whether the command outputs, written by `fn write_nwk_tree()` in the same file, keep it

## Validation

- Parser, writer, Newick round-trip, and Nexus round-trip tests for present, absent, zero, and scientific-notation root lengths
- Root branch lengths round-trip without value loss
- Trees without a root branch length keep their current serialization
- Every `NewickGraph` transformation preserves the field

## Location

- Parser: `fn Builder::finish()` in [packages/util-newick/src/parse.rs](../../packages/util-newick/src/parse.rs)
- Data model: `struct NewickGraph` and `struct NexusTree` in [packages/util-newick/src/types.rs](../../packages/util-newick/src/types.rs)
