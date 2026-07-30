# Point * Function multiplication ignores Function tails

`fn multiply_point_function()` at [packages/treetime-distribution/src/distribution_ops/multiply.rs#L83-L99](../../packages/treetime-distribution/src/distribution_ops/multiply.rs#L83-L99) returns `Empty` when the point falls outside the function's finite grid, even if the function has a `Constant` tail that covers that region.

The product should be `point.amplitude * function.boundary_value` when the point is in the function's tail region, not `Empty`.

Not a production path (inference passes don't produce Point * Function with tails), but inconsistent with the tail-aware contract documented in [kb/decisions/distribution-tails-and-arithmetic.md](../decisions/distribution-tails-and-arithmetic.md).
