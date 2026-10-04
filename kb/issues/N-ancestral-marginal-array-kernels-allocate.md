# Marginal ndarray kernels allocate avoidable intermediates

Dense and sparse marginal inference still create temporaries inside loops that scale with the core inference dimensions: nodes, children, and variable sites.

## Evidence

- Dense backward combination, `fn indexed_node_backward()` [packages/treetime/src/partition/marginal/shared/pass.rs#L106-L108](../../packages/treetime/src/partition/marginal/shared/pass.rs#L106-L108), allocates a `mapv(f64::ln)` array of the full profile for every child before adding it to the running sum.
- Sparse transition statistics, `fn accumulate_site_transition_weighted()` [packages/treetime/src/partition/marginal/sparse/count.rs#L117](../../packages/treetime/src/partition/marginal/sparse/count.rs#L117), allocate a joint `n_states x n_states` array for every variable site of every edge and fill it with nested indexing instead of ndarray broadcasting, `Zip`, and `sum_axis()`.

## Options

- **ndarray operations with borrowed inputs:** use `Zip`, broadcasting, and `sum_axis()`, writing elementwise transforms directly into one result buffer.
- **Reusable scratch arrays with scalar loops:** preallocate buffers while retaining the current indexing. This controls allocation but preserves duplicated indexing and adds mutable scratch-state contracts.

## Recommendation

Express each operation through ndarray and keep a reusable buffer only where measurement shows ndarray still allocates inside the repeated kernel. Any change must keep the summation order, because outputs are compared byte for byte; a reordered sum changes the last bits of the transition counts and of the inferred GTR.

## Required properties

- Whole-array outputs remain equal for finite inputs across binary and multifurcating nodes, dense and sparse partitions, and multiple state counts.
- Existing zero and subnormal finite-value behavior remains unchanged.
- Allocation growth is measured across sites, states, edges, and child count.

## Related issues

- [N-array-owned-signatures-force-projection-copies.md](N-array-owned-signatures-force-projection-copies.md)
