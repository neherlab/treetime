# Tree builder's Fitch refresh materializes full sequences for every node

The Fitch refresh step of the tree builder converts the MAT to a graph, runs `fn create_fitch_partition()` on the full alignment, and reads back a full sequence for every node. The MAT itself stores only edge substitutions.

## Location

`fn fitch_refresh()`: [`packages/treetime/src/build/refresh.rs#L14`](../../packages/treetime/src/build/refresh.rs#L14). It copies each sample sequence, replacing gaps with `N`, and builds the Fitch partition on all sites.

## Impact

Memory and time scale with $n \times L$, not with the number of substitutions:

- On `mpox/clade-ii/2000` (197 kb genomes) the build peaks at 1.2 GB
- The refreshes take 17 of the 21 s of the build
- On `sc2/4500` they take most of the 13.6 s

The refresh is the step with the largest parsimony gain, so it cannot be skipped.

## Fix

Run Fitch only on the sites with a substitution in the MAT or a difference in a query. All other sites are constant across the tree and keep the reference state. This is the parsimony part of [kb/proposals/mutation-first-sequence-representation.md](../proposals/mutation-first-sequence-representation.md). See [kb/proposals/native-parsimony-tree-builder.md](../proposals/native-parsimony-tree-builder.md), axis B5.
