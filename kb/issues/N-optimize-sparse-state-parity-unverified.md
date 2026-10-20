# Sparse branch optimization state parity lacks a coefficient oracle

> [!IMPORTANT]
> **Investigation required.** No failing fixture shows that the stored sparse reference state is stale. Construct the dense-versus-sparse coefficient oracle listed under "Investigation required" and trace any divergence to its root cause.

## Summary

Sparse and dense branch optimization have no coefficient-level oracle proving that their compressed state selection remains equivalent after marginal updates. `VarPos.state` is the exact sparse reference state and cannot be replaced independently with a posterior MAP state.

## Details

`fn get_coefficients()` [packages/treetime/src/partition/optimize/sparse.rs#L11-L84](../../packages/treetime/src/partition/optimize/sparse.rs#L11-L84) uses complete probability profiles for explicit variable sites and compressed profiles plus multiplicities for fixed sites. When one side lacks an explicit variable profile, the state of that side selects its fixed profile: the Fitch substitution states of the edge at that position, or else the stored `VarPos.state` of the other message [packages/treetime/src/partition/optimize/sparse.rs#L29-L55](../../packages/treetime/src/partition/optimize/sparse.rs#L29-L55).

The dense-versus-sparse equivalence tests compare only results, not coefficients: log-likelihoods within 0.5 and branch lengths within 0.05 [packages/treetime/src/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_bounds.rs#L68](../../packages/treetime/src/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_bounds.rs#L68), [packages/treetime/src/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_bounds.rs#L144](../../packages/treetime/src/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_bounds.rs#L144) ([N-optimize-equivalence-bounds-undocumented.md](N-optimize-equivalence-bounds-undocumented.md)).

Changing that key to the current posterior MAP state without rebuilding profiles, reference states, fixed counts, gaps, and unknowns together can select a fixed profile for a different state. No current failing fixture establishes that the stored reference state is stale or that this isolated substitution is valid.

The current evidence does not establish that the reference state should contain a posterior MAP state or that any stored representation component is stale.

## Investigation required

- Construct a dense-versus-sparse oracle for every coefficient and optimized edge length after repeated marginal updates.
- Include explicit variable profiles, fixed rows, gaps, unknowns, and changed MAP assignments.
- Trace the first divergence to the representation component that is stale.
- If stale compression is confirmed, specify how all coupled compressed state is rebuilt atomically.

Implementation waits until the coefficient oracle reproduces a defect and identifies its root cause.
