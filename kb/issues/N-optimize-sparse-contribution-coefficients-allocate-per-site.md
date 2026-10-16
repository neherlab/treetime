# Sparse branch-length coefficients allocate per site

`fn get_coefficients()` [packages/treetime/src/partition/optimize/sparse.rs#L11-L84](../../packages/treetime/src/partition/optimize/sparse.rs#L11-L84) builds the branch-length contribution of one edge. For every explicit position and every fixed state it computes `aview1(child).dot(&gtr.v) * aview1(parent).dot(&gtr.v_inv.t())` and stores the result as the `Array1` coefficients of a `SiteContribution`, so each position allocates the two `dot` results on the heap.

On `rsv/a/500`, sparse `timetree` with one thread, these allocations are 0.82 M of the 6.33 M allocations of the whole run (13 %), measured with the heap-profiling build on `rust` at `d88dc18f`.

## Fix direction

Keep the coefficients inline, as the sparse marginal messages keep their profiles (`type SiteProfile` in [packages/treetime/src/partition/storage/sparse.rs](../../packages/treetime/src/partition/storage/sparse.rs)), and write both products into inline buffers with `ndarray::linalg::general_mat_vec_mul`. Keep the evaluation order of each coefficient (one BLAS matrix-vector product per factor, then the elementwise product), so branch lengths stay byte-identical.

## Related issues

- [N-optimize-edge-contribution-copies-gtr.md](N-optimize-edge-contribution-copies-gtr.md): the same function copies the GTR model for every edge
- [N-ancestral-marginal-array-kernels-allocate.md](N-ancestral-marginal-array-kernels-allocate.md): per-site temporaries of the marginal kernels
