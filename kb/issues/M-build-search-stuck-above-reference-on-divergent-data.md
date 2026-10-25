# Built trees stay 0.7% above the reference tree on divergent data

On `lassa/L/500` the best built tree (neighbor-joining start, SPR, rebuilds of radius 3 or 5) scores 75,674 to 75,989 Fitch substitutions and misses 82 to 95 of the 411 reference splits. The reference tree, refined by the same moves, scores 75,258. That is lower than every built tree, and it keeps 406 of its splits. So the objective is not the limit: the search ends in worse local optima, and 66 of the missed splits sit on reference branches longer than 10 substitutions.

Building from slow sites first, or giving them more weight, makes the trees worse. Half weight for the fastest sites, and a parsimony ratchet (random reweighting of 25% of sites, eight iterations), close 0.1 points ([kb/reports/feat-tree-infer-prototype.md](../reports/feat-tree-infer-prototype.md), section "Site partitions and reweighting on divergent data").

## Location

`fn refine_mat()`: [`packages/treetime/src/build/pipeline.rs`](../../packages/treetime/src/build/pipeline.rs)

## Impact

Trees for divergent pathogens (Lassa and other arenaviruses, HIV-like diversity) differ from ML trees in deep splits, which also changes the root and the clock estimates.

## Fix

Add site weights to the MAT score, run the ratchet inside the pipeline instead of through rewritten alignments, and measure the score against the number of iterations with stronger perturbations. Compare with several neighbor-joining starts on resampled sites. Accept when `lassa/L/500` ends within 0.2% of the refined reference tree.
