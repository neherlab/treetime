# Placement bound counts every substitution below a node

The branch-and-bound placement prunes the subtree of $v$ when $c(v) - d(v)$ is not below the best cost, where $d(v)$ is the longest substitution path below $v$. A query can gain only at sites where it differs from the state of $v$. Substitutions at other sites can only raise its cost, yet they still loosen the bound.

## Location

`fn Mat.depth_below()` and `fn Mat.raise_depth()`: [`packages/treetime/src/build/mat.rs#L84-L105`](../../packages/treetime/src/build/mat.rs#L84-L105), used by the search in [`packages/treetime/src/build/placement.rs`](../../packages/treetime/src/build/placement.rs)

## Impact

The search visits the same set of nodes from any starting point. Starting from reference nodes saved only 7 to 30% of the visits ([kb/reports/feat-tree-infer-prototype.md](../reports/feat-tree-infer-prototype.md)). The visit count grows with the tree and is expected to dominate at 10^5 samples. There is no measurement at that scale yet.

## Fix

Bound by the number of the query's differing sites that are substituted below $v$, for example with a per-node set or sketch of the positions substituted below it. The bound must stay valid after insertions and detachments. The existing oracle tests (branch-and-bound equal to the full scan) cover correctness.
