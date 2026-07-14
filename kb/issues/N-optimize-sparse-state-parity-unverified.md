# Sparse branch optimization state parity lacks a coefficient oracle

## Summary

Sparse and dense branch optimization have no coefficient-level oracle proving that their compressed state selection remains equivalent after marginal updates. `VarPos.state` is the exact sparse reference state and cannot be replaced independently with a posterior MAP state.

## Details

`fn get_coefficients()` [`packages/treetime/src/partition/optimize_sparse.rs#L45-L118`](../../packages/treetime/src/partition/optimize_sparse.rs#L45-L118) uses complete probability profiles for explicit variable sites and compressed profiles plus multiplicities for fixed sites. The reference state selects the fixed profile when one side lacks an explicit variable profile.

Changing that key to the current posterior MAP state without rebuilding profiles, reference states, `fixed_counts`, gaps, and unknowns together can select a fixed profile for a different state. No current failing fixture establishes that the stored reference state is stale or that this isolated substitution is valid.

The current evidence does not establish that the reference state should contain a posterior MAP state or that any stored representation component is stale.

## Investigation required

- Construct a dense-versus-sparse oracle for every coefficient and optimized edge length after repeated marginal updates.
- Include explicit variable profiles, fixed rows, gaps, unknowns, and changed MAP assignments.
- Trace the first divergence to the representation component that is stale.
- If stale compression is confirmed, specify how all coupled compressed state is rebuilt atomically.

No implementation ticket is ready until the coefficient oracle reproduces a defect and identifies its root cause.
