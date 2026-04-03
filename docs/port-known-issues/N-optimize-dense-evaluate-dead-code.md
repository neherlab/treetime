# Dead `optimize_dense::evaluate()` function

`optimize_dense::evaluate()` at [packages/treetime/src/commands/optimize/optimize_dense.rs#L39-L57](../../packages/treetime/src/commands/optimize/optimize_dense.rs#L39-L57) duplicates `evaluate_dense_contribution_impl()` from `optimize_dense_eval.rs` with a slightly different structure (always computes derivatives, no option to skip). It is not called from production code - only from two test files:

- `test_coefficient_extraction_dense_basic.rs`
- `test_coefficient_extraction_dense_properties.rs`

## Impact

Negligible. Dead production code that increases maintenance surface. Tests should use `evaluate_dense_contribution()` from `optimize_dense_eval.rs` instead.

## Fix

Remove `optimize_dense::evaluate()`, update the two test files to import `evaluate_dense_contribution` from `optimize_dense_eval`.
