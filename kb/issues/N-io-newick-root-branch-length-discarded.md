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

1. Add `root_branch_length: Option<f64>` to `NewickGraph`
2. Store the parsed root branch length in `fn visit_root_branch()` instead of discarding it
3. Emit the stored value after the root subtree in the Newick writer `fn newick_to_writer()` [`packages/util-newick/src/write.rs#L7`](../../packages/util-newick/src/write.rs#L7)
4. Carry the field through `NexusTree` records and every conversion that constructs or consumes `NewickGraph`

## Validation

- Parser, writer, Newick round-trip, and Nexus round-trip tests for present, absent, zero, and scientific-notation root lengths
- Root branch lengths round-trip without value loss
- Trees without a root branch length keep their current serialization
- Every `NewickGraph` transformation preserves the field

## Location

- Parser: `fn visit_root_branch()` [`packages/util-newick/src/parse.rs#L77-L105`](../../packages/util-newick/src/parse.rs#L77-L105)
- Data model: `struct NewickGraph` [`packages/util-newick/src/types.rs#L30-L35`](../../packages/util-newick/src/types.rs#L30-L35), `struct NexusTree` [`packages/util-newick/src/types.rs#L24-L27`](../../packages/util-newick/src/types.rs#L24-L27)
