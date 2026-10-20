# Manual index loops replace available ndarray operations

Several numerical routines build or reduce arrays with nested index loops where an ndarray operation or a `treetime-utils` array helper expresses the same computation. The loops hide the array semantics (outer products, axis sums, reversal), repeat bounds checks per element, and diverge from the project rule to stay inside `ndarray` for computation.

## Locations

- `fn get_branch_mutation_matrix()` [packages/treetime/src/gtr/infer_gtr/common.rs#L130-L159](../../packages/treetime/src/gtr/infer_gtr/common.rs#L130-L159): fills a sites x states x states joint array with a triple loop and normalizes each site with a second double loop, instead of broadcasting `msg_to_parent`, `exp_qt`, and `msg_to_child` and dividing by the per-site `sum_axis()`
- `fn accumulate_mutation_counts()` [packages/treetime/src/gtr/infer_gtr/common.rs#L161-L192](../../packages/treetime/src/gtr/infer_gtr/common.rs#L161-L192): triple-nested loops that compute `nij` and `Ti` as sums over sites and states, replaceable by `sum_axis()`
- `fn jtt92()` [packages/treetime/src/gtr/get_gtr.rs#L355-L363](../../packages/treetime/src/gtr/get_gtr.rs#L355-L363): builds `W` from `Q` and `pi` with a double loop instead of an elementwise expression over `Q` and the outer ratio of `pi`
- `fn interp_nonuniform()` [packages/treetime-grid/src/interp_nonuniform.rs#L26-L41](../../packages/treetime-grid/src/interp_nonuniform.rs#L26-L41): allocates zeros and assigns each element in a loop instead of `Array1::from_shape_fn`
- `fn negate_arg_inplace()` [packages/treetime-grid/src/grid_fn.rs#L364-L367](../../packages/treetime-grid/src/grid_fn.rs#L364-L367): manual swap loop instead of `reverse_inplace()` from `treetime-utils` [packages/treetime-utils/src/array/ndarray.rs#L113](../../packages/treetime-utils/src/array/ndarray.rs#L113)
- `fn compute_relative_error_statistics()` [packages/treetime-validation/src/testing/metrics/aggregate/domain_agreement/error_stats.rs#L31-L39](../../packages/treetime-validation/src/testing/metrics/aggregate/domain_agreement/error_stats.rs#L31-L39): copies filtered ndarray elements into `Vec`s for the statistics instead of computing them over ndarray

## Fix direction

Express each computation through ndarray (`Zip`, broadcasting, `sum_axis()`, `from_shape_fn`) or the matching `treetime-utils` helper. Outputs are compared byte for byte, so keep the summation order of each reduction or confirm that the reordered sum leaves the golden-master and smoke outputs unchanged.

## Related

- [N-ancestral-sparse-transition-count-allocates-per-site.md](N-ancestral-sparse-transition-count-allocates-per-site.md): the sparse counterpart of `get_branch_mutation_matrix()` uses the same per-site joint-array loop
