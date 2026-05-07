# Dense/sparse equivalence test bounds undocumented

The equivalence tests in `test_dense_sparse_equivalence_bounds.rs` use integration-level behavioral bounds (log-LH difference < 0.5, total tree length difference < 0.1, per-edge branch length difference < 0.05) without inline documentation of their derivation or expected tightening criteria.

## Impact

The bounds are appropriate for comparing two different algorithmic representations converging to similar ML estimates. They are not floating-point precision tolerances. The lack of documentation makes it harder to evaluate whether a regression changed the bounds or the bounds were always loose.

## Proposed solution

Add inline comments explaining why each bound is set to its current value and under what conditions it could be tightened.

## Related issues

- Source: [N-optimize-equivalence-bounds-undocumented.md](../issues/N-optimize-equivalence-bounds-undocumented.md) -- delete after full resolution
