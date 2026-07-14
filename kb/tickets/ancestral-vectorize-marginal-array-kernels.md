# Vectorize marginal ndarray kernels

Replace deep clones and full-size intermediates in dense and sparse marginal passes with borrowed, shape-aware ndarray operations.

## Required changes

- Accumulate dense child logarithms and cavity division with `Zip` in [packages/treetime/src/partition/marginal_core.rs#L151](../../packages/treetime/src/partition/marginal_core.rs#L151) and [packages/treetime/src/partition/marginal_core.rs#L244](../../packages/treetime/src/partition/marginal_core.rs#L244).
- Accept iterators of borrowed sparse messages instead of cloning profile maps in [packages/treetime/src/partition/marginal_passes.rs#L164](../../packages/treetime/src/partition/marginal_passes.rs#L164) and [packages/treetime/src/partition/marginal_passes.rs#L327](../../packages/treetime/src/partition/marginal_passes.rs#L327).
- Use broadcasting or `Zip` for joint transition matrices and `sum_axis` for marginals in [packages/treetime/src/partition/marginal_sparse.rs#L221](../../packages/treetime/src/partition/marginal_sparse.rs#L221).
- Preserve deterministic ordering and existing finite-input semantics.

## Validation

- Whole-array dense and sparse equivalence for binary and multifurcating nodes.
- Zero, subnormal, and representative positive finite inputs with exact equality where operation order is unchanged and `max_ulps = 10` where ndarray changes only floating-point evaluation order.
- Allocation benchmarks scaling by sites, states, edges, and children.
- Full lint and test suite.

## Related issues

- Source: [kb/issues/N-ancestral-marginal-array-kernels-allocate.md](../issues/N-ancestral-marginal-array-kernels-allocate.md)
