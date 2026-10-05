# Dense/sparse equivalence test bounds undocumented

> [!WARNING]
> **Needs review.** The project bans code comments in Rust source, test files included, so the bounds cannot be documented with inline comments. The documentation needs another carrier, for example named constants whose names state the meaning of each bound, or a KB entry that records the derivation and the tightening criteria.

The equivalence tests in [`packages/treetime/src/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_bounds.rs`](../../packages/treetime/src/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_bounds.rs) use integration-level behavioral bounds without documentation of their derivation or expected tightening criteria:

- log-LH difference < 0.5 [test_dense_sparse_equivalence_bounds.rs#L71](../../packages/treetime/src/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_bounds.rs#L71)
- total tree length difference < 0.1 [test_dense_sparse_equivalence_bounds.rs#L141](../../packages/treetime/src/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_bounds.rs#L141)
- per-edge branch length difference < 0.05 [test_dense_sparse_equivalence_bounds.rs#L151](../../packages/treetime/src/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_bounds.rs#L151)

## Impact

The bounds are appropriate for comparing two different algorithmic representations converging to similar ML estimates. They are not floating-point precision tolerances. The lack of documentation makes it harder to evaluate whether a regression changed the bounds or the bounds were always loose.

## Proposed solution

Record why each bound has its current value and under what conditions it could be tightened, in a form the comment ban permits.
