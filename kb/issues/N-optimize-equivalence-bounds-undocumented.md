# Dense/sparse equivalence test bounds undocumented

The equivalence tests in [`packages/treetime/src/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_bounds.rs`](../../packages/treetime/src/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_bounds.rs) use integration-level behavioral bounds without documentation of their derivation or expected tightening criteria:

- log-LH difference < 0.5 [test_dense_sparse_equivalence_bounds.rs#L68](../../packages/treetime/src/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_bounds.rs#L68)
- total tree length difference < 0.1 [test_dense_sparse_equivalence_bounds.rs#L134](../../packages/treetime/src/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_bounds.rs#L134)
- per-edge branch length difference < 0.05 [test_dense_sparse_equivalence_bounds.rs#L144](../../packages/treetime/src/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_bounds.rs#L144)

## Impact

The bounds are appropriate for comparing two different algorithmic representations converging to similar ML estimates. They are not floating-point precision tolerances. The lack of documentation makes it harder to evaluate whether a regression changed the bounds or the bounds were always loose.

## Proposed solution

Record why each bound has its current value and under what conditions it could be tightened. The project bans code comments in Rust source, test files included, so the carrier is a named constant whose name states the meaning of each bound, the assertion message, or a KB entry with the derivation and the tightening criteria.

The bounds also conflict with the project rule that a numerical error above 1e-6 in a test assertion is a bug, unless the dense and sparse paths differ by design for the tested model (the sparse path approximates per-site rates with a scalar rate at fixed sites).

## Proposed implementation plan

1. Measure: print to stderr the three differences (log-LH, total length, maximum per-edge length) of each test in `test_dense_sparse_equivalence_bounds.rs` on the current code, and record which inputs use per-site rates, gamma categories, or indels.
2. Classify each test: when the model has no approximation in sparse mode (no per-site rates), dense and sparse are independent oracles of one likelihood, and the bound becomes the project tolerance (`assert_abs_diff_eq!` or `assert_ulps_eq!` with the tightest passing `1e-N`, at most `1e-6`). A measured difference above that is a defect: open a new issue with the reproduction instead of keeping the bound.
3. For a test whose model has a documented sparse approximation, derive the bound from that approximation (cite `kb/decisions/` or the project rule that approves it) and state the derivation in the assertion message, or in a constant whose name states the meaning, because the project bans comments in Rust source.
4. Replace the `assert!(diff < x)` checks with `approx` assertions on whole values (the branch-length arrays compared as arrays), so a failure prints both arrays.
5. Validate with `just test-rs -E 'package(treetime)' test_dense_sparse_equivalence`, then see red: set one expected tolerance to `1e-20` and confirm the failure message, then revert.
