# Dense and sparse branch-length evaluators duplicate the site iteration

`fn evaluate_dense_contribution()` and `fn evaluate_sparse_contribution()` each build a `(multiplicity, coefficients)` iterator and pass it to `fn evaluate_site_contributions()` [packages/treetime/src/optimize/dense_eval.rs#L6](../../packages/treetime/src/optimize/dense_eval.rs#L6), [packages/treetime/src/optimize/sparse_eval.rs#L6](../../packages/treetime/src/optimize/sparse_eval.rs#L6). `OptimizationContribution::sites()` already builds the same iterator for both variants [packages/treetime/src/partition/optimize/contribution.rs#L30](../../packages/treetime/src/partition/optimize/contribution.rs#L30), so the two files repeat logic that the enum owns.

## Proposed solution

Let `fn evaluate_mixed_impl()` [packages/treetime/src/optimize/likelihood.rs](../../packages/treetime/src/optimize/likelihood.rs) call `evaluate_site_contributions(contribution.sites(), contribution.rate_eigvals(), ...)` directly, and delete `dense_eval.rs` and `sparse_eval.rs`. The tests that call the two wrappers (the `test_coefficient_extraction_*` modules) then wrap their contributions in `OptimizationContribution` and evaluate through the same entry point.
