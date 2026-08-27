# Distribution operation error-message casing mismatch fails divide/convolve tests

## Symptom

On `dev`, the `treetime-distribution` test suite is red. The parameterized cases in `packages/treetime-distribution/src/distribution_ops/__tests__/test_divide.rs` (`test_divide_formula_combinations_return_errors`) and the matching `test_convolve.rs` cases (`test_convolve_formula_combinations_return_errors`) fail because the produced error message and the expected string disagree on the casing of the operand kind names.

Example (`case_1_formula_empty`):

- expected: `Cannot divide Formula by Empty: operation not implemented`
- actual: `Cannot divide formula by empty: operation not implemented`

All nine divide cases and the corresponding convolve cases fail the same way.

## Root cause

`889a8790 refactor(distribution): simplify operation variant dispatch` changed how the unsupported-operation error renders the operand kind, so the message now lowercases the kind name (`formula`, `empty`) while the tests still expect the capitalized variant names (`Formula`, `Empty`). Either the error text or the test expectations were updated without the other.

## Scope

Presentation/test defect: the mismatch is in a user-facing error string, with no demonstrated effect on numerical results. It does block the `treetime-distribution` test suite.

## Fix approach

Decide the intended rendering of operand kind names in distribution operation errors (capitalized variant names or lowercase words), then align the `Display`/error construction and the divide/convolve test expectations to match. Verify no other distribution error messages depend on the old casing.
