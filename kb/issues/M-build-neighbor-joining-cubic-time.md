# Neighbor-joining start of the tree builder takes cubic time

The tree builder's distance start runs neighbor joining with a full scan of the Q-criterion in every join. It takes $O(n^3)$ time and $O(n^2)$ memory. It was not run on `sc2/4500` (5100 samples) for this reason. On `lassa/L/500` it is the start that matters: it lowers the excess over the reference tree from +3.1% to +0.96% after the same refinement.

## Location

- `fn neighbor_joining()`: [`packages/treetime/src/build/nj.rs#L31-L83`](../../packages/treetime/src/build/nj.rs#L31-L83). It recomputes the row sums and scans all pairs in each of the $n - 3$ joins
- `fn jc_distances()`: [`packages/treetime/src/build/nj.rs#L12-L29`](../../packages/treetime/src/build/nj.rs#L12-L29). It fills the full $n \times n$ matrix

## Impact

Divergent datasets with more than a few thousand samples cannot use the distance start.

## Fix

Options from [kb/proposals/native-parsimony-tree-builder.md](../proposals/native-parsimony-tree-builder.md), axis B2:

- Rapid NJ with bounds on the Q-criterion: the same tree, with most pairs skipped
- NJ on a random subsample of a few hundred samples, then placement of the remaining samples. Untested; whether a subsample keeps the deep structure needs a check on `lassa/L/500` and on rate-variation simulations

Any faster method needs a test against the current implementation on an additive metric, where NJ recovers the true tree.
