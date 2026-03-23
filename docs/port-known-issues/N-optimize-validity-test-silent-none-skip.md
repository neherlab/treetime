# Dense/sparse equivalence validity tests silently skip None branch lengths

The validity tests `test_dense_optimization_produces_valid_results` and `test_sparse_optimization_produces_valid_results` use `if let Some(bl)` to guard branch-length assertions. Since `run_optimize_mixed()` sets `Some(...)` on every edge, `None` cannot occur post-optimization. If a regression caused `None`, the test would pass silently instead of failing.

## Location

[`packages/treetime/src/commands/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_validity.rs#L47`](../../packages/treetime/src/commands/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_validity.rs#L47) and line 88.

## Proposed solution

Replace `if let Some(bl)` with `let bl = bl.expect("branch length must be set after optimization")` or `.unwrap()`.
