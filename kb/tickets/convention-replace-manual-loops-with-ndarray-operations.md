# Replace 7 manual loops with ndarray operations

Seven locations use manual loops where ndarray operations would be clearer and more idiomatic.

## Instances

1. `packages/treetime/src/gtr/infer_gtr/common.rs#L100-L129` -- triple-nested loop replaceable by `sum_axis()`
2. `packages/treetime/src/gtr/infer_gtr/common.rs#L59-L79` -- manual outer product
3. `packages/treetime/src/gtr/get_gtr.rs#L398-L409` -- manual W computation
4. `packages/treetime-grid/src/interp_nonuniform.rs#L38-L53` -- zeros + loop instead of `from_shape_fn`
5. `packages/treetime-grid/src/grid_fn.rs#L286-L289` -- manual swap loop instead of `reverse_inplace()`
6. `packages/treetime-validation/src/testing/metrics/aggregate/domain_agreement/error_stats.rs#L43-L44` -- ndarray to Vec conversion

## Related issues

- Source: [kb/issues/N-code-quality-conventions.md](../issues/N-code-quality-conventions.md) -- delete after full resolution
