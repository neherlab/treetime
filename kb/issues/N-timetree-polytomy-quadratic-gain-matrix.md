# Polytomy resolution recomputes full gain matrix each iteration

`resolve_single_polytomy` recomputes the full O(n^2) pairwise gain matrix from scratch after each merge by re-reading the graph topology. v0 updates the matrix incrementally: after merging a pair, it removes their rows/cols and computes only O(n) new entries for the merged node vs remaining children.

v0: `packages/legacy/treetime/treetime/treetime.py:838-856`
v1: `packages/treetime/src/timetree/optimization/polytomy.rs`

For a k-way polytomy, current cost is sum(n^2) for n from k down to 3. Incremental cost is k\*(k-1)/2 initial + sum(n-2) for subsequent merges.
