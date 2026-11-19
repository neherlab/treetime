# Manual index loops replace available ndarray operations

Several numerical routines build or reduce arrays with nested index loops where an ndarray operation or a `treetime-utils` array helper expresses the same computation. The loops hide the array semantics (outer products, axis sums, reversal), repeat bounds checks per element, and diverge from the project rule to stay inside `ndarray` for computation.

## Locations

- `fn get_branch_mutation_matrix()` [packages/treetime/src/gtr/infer_gtr/common.rs#L250-L279](../../packages/treetime/src/gtr/infer_gtr/common.rs#L250-L279): fills a sites x states x states joint array with a triple loop and normalizes each site with a second double loop, instead of broadcasting `msg_to_parent`, `exp_qt`, and `msg_to_child` and dividing by the per-site `sum_axis()`
- `fn accumulate_mutation_counts()` [packages/treetime/src/gtr/infer_gtr/common.rs#L281-L312](../../packages/treetime/src/gtr/infer_gtr/common.rs#L281-L312): triple-nested loops that compute `nij` and `Ti` as sums over sites and states, replaceable by `sum_axis()`
- `fn jtt92()` [packages/treetime/src/gtr/get_gtr.rs#L544-L551](../../packages/treetime/src/gtr/get_gtr.rs#L544-L551): builds `W` from `Q` and `pi` with a double loop instead of an elementwise expression over `Q` and the outer ratio of `pi`
- `fn negate_arg_inplace()` [packages/treetime-grid/src/grid_fn.rs#L361-L368](../../packages/treetime-grid/src/grid_fn.rs#L361-L368): manual swap loop instead of `reverse_inplace()` from `treetime-utils` [packages/treetime-utils/src/array/ndarray.rs#L132](../../packages/treetime-utils/src/array/ndarray.rs#L132)
- `fn compute_relative_error_statistics()` [packages/treetime-validation/src/testing/metrics/aggregate/domain_agreement/error_stats.rs#L31-L41](../../packages/treetime-validation/src/testing/metrics/aggregate/domain_agreement/error_stats.rs#L31-L41): copies filtered ndarray elements into `Vec`s for the statistics instead of computing them over ndarray

## Fix direction

Express each computation through ndarray (`Zip`, broadcasting, `sum_axis()`, `from_shape_fn`) or the matching `treetime-utils` helper. Outputs are compared byte for byte, so keep the summation order of each reduction or confirm that the reordered sum leaves the golden-master and smoke outputs unchanged.

## Proposed implementation plan

All five locations are current (`fn jtt92()` starts at `get_gtr.rs` L507; its loop is at L543-L551).

1. `fn negate_arg_inplace()`: replace the swap loop with `reverse_inplace(&mut self.y)`
2. `fn jtt92()`: build `W` with `Array2::from_shape_fn((n, n), |(i, j)| ...)`, zero on the diagonal
3. `fn get_branch_mutation_matrix()`: broadcast `msg_to_parent`, `exp_qt`, and `msg_to_child` into the sites x states x states array, and divide by the per-site `sum_axis()` with an inserted axis
4. `fn accumulate_mutation_counts()`: compute `nij` and `Ti` with `sum_axis(Axis(0))` over the joint array
5. `fn compute_relative_error_statistics()`: compute the statistics over the filtered ndarray view or a `Zip` instead of `Vec` copies
6. Tests: the existing tests of these functions stay unchanged and must pass (oracle: their current expected values, which come from v0 for the GTR inference); for a reordered reduction in steps 3 and 4, compare with `assert_ulps_eq!` and confirm that the golden-master tests and `just smoke` on `data/flu/h3n2/20` give unchanged outputs
7. Validation: `just check`, `just test-rs`, and `just smoke` for steps 3 and 4
