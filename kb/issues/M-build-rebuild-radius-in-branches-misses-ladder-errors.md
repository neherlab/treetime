# Local rebuilds measure their radius in branches, which misses errors along ladders

The tree builder's local rebuild cuts out the ball of radius $r$ branches around each internal node. On `sc2/4500`, 41% of the misplaced tips of conflicting splits lie more than 10 branches from the split, but only 11 substitutions away (median). They lie along chains of branches with one substitution each, so no ball of 3 to 5 branches contains them. Rebuilds removed 3 to 10% of the conflicting splits on `sc2/4500`, against 38 to 52% on `rsv/a/2000` and `dengue/2000` ([kb/reports/feat-tree-infer-prototype.md](../reports/feat-tree-infer-prototype.md)).

At the same time, large polytomies push about a third of the balls of radius 5 on `sc2/4500` and `mpox/clade-ii/2000` over the size limit of 48 pseudo-leaves.

## Location

`fn Ball::around()`: [`packages/treetime/src/build/rebuild.rs#L38-L71`](../../packages/treetime/src/build/rebuild.rs#L38-L71)

## Impact

`sc2/4500` stays +0.29% above the reference tree, the largest gap among the low-divergence datasets.

## Fix

Define the ball by substitution distance, capped by the number of pseudo-leaves, and compare it with the branch radius on `sc2/4500` and the simulations. See [kb/proposals/native-parsimony-tree-builder.md](../proposals/native-parsimony-tree-builder.md), axis B4.
