# Neighbor-joining start of the tree builder needs the full distance matrix

The tree builder's distance start computes and stores all pairwise Jukes-Cantor distances before it joins. Time and memory are $O(n^2)$: about 0.8 GB at $10^4$ samples and 80 GB at $10^5$.

The join search itself is no longer cubic. With top hits ($m \approx \sqrt n$, FastTree's heuristic), it costs $O(n^3/m)$ and gives trees identical to exact neighbor joining on `lassa/L/500` and `dengue/2000`. `sc2/4500` takes 17 s ([kb/reports/feat-tree-infer-prototype.md](../reports/feat-tree-infer-prototype.md#fast-and-local-neighbor-joining)).

## Location

- `fn jc_distances()`: [`packages/treetime/src/build/nj.rs#L12-L29`](../../packages/treetime/src/build/nj.rs#L12-L29), which fills the full $n \times n$ matrix
- `fn neighbor_joining()` in the same file, which updates that matrix in place

## Impact

The neighbor-joining start, the better start on `lassa/L/500` and `sc2/4500`, cannot run on more than about $10^4$ samples.

## Fix

Options from [kb/proposals/native-parsimony-tree-builder.md](../proposals/native-parsimony-tree-builder.md), axis B2:

- Profiles of joined nodes instead of a matrix, as in FastTree: $O(n \sqrt n)$ memory
- Neighbor joining on a subsample, then placement of the remaining samples. Untested

Any replacement needs a test against the current implementation on an additive metric, where neighbor joining recovers the true tree.
