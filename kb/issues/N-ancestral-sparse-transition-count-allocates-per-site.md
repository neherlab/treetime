# Sparse marginal transition counting allocates a joint array per site

Sparse transition statistics create a temporary inside a loop that scales with the core inference dimensions: edges and variable sites.

## Evidence

- `fn accumulate_site_transition_weighted()` [packages/treetime/src/partition/marginal/sparse/count.rs#L91-L132](../../packages/treetime/src/partition/marginal/sparse/count.rs#L91-L132) allocates a joint `n_states x n_states` array for every variable site of every edge and fills, normalizes, and sums it with nested index loops instead of ndarray broadcasting, `Zip`, and `sum_axis()`.

The dense backward combination does not have this problem: `fn indexed_node_backward()` [packages/treetime/src/partition/marginal/shared/pass.rs#L162-L165](../../packages/treetime/src/partition/marginal/shared/pass.rs#L162-L165) allocates one log array per node and adds the other children into it in place.

## Options

- **ndarray operations with borrowed inputs:** use `Zip`, broadcasting, and `sum_axis()`, writing elementwise transforms directly into one result buffer.
- **Reusable scratch arrays with scalar loops:** preallocate buffers while retaining the current indexing. This controls allocation but preserves duplicated indexing and adds mutable scratch-state contracts.

## Recommendation

Express the operation through ndarray and keep a reusable buffer only where measurement shows ndarray still allocates inside the repeated kernel. Any change must keep the summation order, because outputs are compared byte for byte; a reordered sum changes the last bits of the transition counts and of the inferred GTR.

## Required properties

- Transition counts remain equal for finite inputs across multiple state counts.
- Existing zero and subnormal finite-value behavior remains unchanged.
- Allocation growth is measured across sites, states, and edges.

## Related issues

- [N-array-owned-signatures-force-projection-copies.md](N-array-owned-signatures-force-projection-copies.md)
- [N-ancestral-sparse-profiles-allocate-for-large-alphabets.md](N-ancestral-sparse-profiles-allocate-for-large-alphabets.md): per-site profile allocations of sparse messages
