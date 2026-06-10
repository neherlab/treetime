# Root branch length silently discarded during Newick parsing

Severity: N (negligible)
Category: I/O
Crate: `util-newick`

## Description

The Newick parser accepts root-level branch lengths (`(A:0.1,B:0.2):0.5;`) but silently discards the root's `:0.5`. `NewickGraph` has no field to represent root branch lengths since the root node has no incoming edge.

Round-trip of trees with root branch lengths loses data: `(A:0.1,B:0.2):0.5;` parses and writes back as `(A:0.1,B:0.2);`.

## Impact

Root branch lengths appear in some tool outputs (IQ-TREE, BEAST) and can carry evolutionary distance information (e.g., stem branch above displayed root). The loss is silent with no warning to callers.

## Fix

Add `root_branch_length: Option<f64>` to `NewickGraph`. Store the parsed root branch length there instead of discarding it. Update the writer to emit it after the root subtree.

## Location

- Parser: `packages/util-newick/src/parse.rs` `visit_root_branch()` line ~93
- Data model: `packages/util-newick/src/types.rs` `NewickGraph`
