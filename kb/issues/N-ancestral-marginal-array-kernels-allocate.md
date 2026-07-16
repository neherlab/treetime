# Marginal ndarray kernels allocate avoidable intermediates

Dense and sparse marginal inference repeatedly creates full-profile temporaries and deep-clones messages for read-only combination. These allocations occur inside node, edge, child, or variable-site loops and scale with the core inference dimensions.

## Evidence

- Dense backward combination [packages/treetime/src/partition/marginal_core.rs#L151](../../packages/treetime/src/partition/marginal_core.rs#L151) allocates a `mapv(f64::ln)` array for every additional child.
- Indexed sparse passes [packages/treetime/src/partition/marginal_passes.rs#L164](../../packages/treetime/src/partition/marginal_passes.rs#L164) and [packages/treetime/src/partition/marginal_passes.rs#L327](../../packages/treetime/src/partition/marginal_passes.rs#L327) clone complete `SparseSeqDistribution` messages even though combination only reads them.
- Dense and sparse cavity calculations [packages/treetime/src/partition/marginal_core.rs#L244](../../packages/treetime/src/partition/marginal_core.rs#L244) and [packages/treetime/src/partition/marginal_passes.rs#L228](../../packages/treetime/src/partition/marginal_passes.rs#L228) allocate a sanitized denominator and then a separate quotient.
- Sparse transition statistics [packages/treetime/src/partition/marginal_sparse.rs#L221](../../packages/treetime/src/partition/marginal_sparse.rs#L221) manually construct joint values and marginals with nested indexing instead of ndarray broadcasting, `Zip`, and `sum_axis()`.

## Options

- **ndarray operations with borrowed inputs:** use `Zip`, broadcasting, `sum_axis()`, and iterators of borrowed messages. Array shapes remain checked by ndarray and elementwise transforms can write directly to one result buffer.
- **Reusable scratch arrays with scalar loops:** preallocate buffers while retaining the current indexing. This controls allocation but preserves duplicated indexing and adds mutable scratch-state contracts.

## Recommendation

Express each mathematical operation through ndarray and borrow message inputs; retain a reusable buffer only where measurement shows ndarray still allocates inside the repeated kernel.

## Required properties

- Whole-array outputs remain equal for finite inputs across binary and multifurcating nodes, dense and sparse partitions, and multiple state counts.
- Existing zero and subnormal finite-value behavior remains unchanged.
- Sparse message ordering remains deterministic.
- Allocation growth is measured across sites, states, edges, and child count.

## Related issues

- [N-array-owned-signatures-force-projection-copies.md](N-array-owned-signatures-force-projection-copies.md)
