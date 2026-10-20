# Dense and sparse branch-length evaluators duplicate the site iteration

`fn evaluate_dense_contribution()` and `fn evaluate_sparse_contribution()` each build a `(multiplicity, coefficients)` iterator and pass it to `fn evaluate_site_contributions()` [packages/treetime/src/optimize/dense_eval.rs#L6](../../packages/treetime/src/optimize/dense_eval.rs#L6), [packages/treetime/src/optimize/sparse_eval.rs#L6](../../packages/treetime/src/optimize/sparse_eval.rs#L6). The private `OptimizationContribution::sites()` already builds the same iterator for both variants [packages/treetime/src/partition/optimize/contribution.rs#L37-L49](../../packages/treetime/src/partition/optimize/contribution.rs#L37-L49), so the two files repeat logic that the enum owns. `fn sum_partition_metrics()` [packages/treetime/src/optimize/likelihood.rs#L42-L56](../../packages/treetime/src/optimize/likelihood.rs#L42-L56) dispatches to the two wrappers.

## Proposed solution

Make `sites()` and `rate_eigvals()` `pub(crate)`, let `fn sum_partition_metrics()` call `evaluate_site_contributions(contribution.sites(), contribution.rate_eigvals(), ...)` directly, and delete `dense_eval.rs` and `sparse_eval.rs`. The tests that call the two wrappers (the `test_coefficient_extraction_*` modules) then wrap their contributions in `OptimizationContribution` and evaluate through the same entry point.
