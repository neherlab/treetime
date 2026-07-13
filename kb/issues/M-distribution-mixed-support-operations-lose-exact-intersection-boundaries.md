# Mixed-support distribution operations lose exact intersection boundaries

Several distribution multiplication and division paths do not preserve the analytical intersection of operand supports. A range boundary that falls between function knots is snapped to an existing knot or ignored, changing the domain on which the result is defined.

## Affected operations

### Range multiplied by Function

[`multiply_range_function()`](../../packages/treetime-distribution/src/distribution_ops/multiply.rs) filters the function's existing knots using the range bounds, then sets the result bounds to the first and last retained knots. For a range such as `[1.5, 3.5]` and a function with integer knots, the result becomes `[2.0, 3.0]` instead of the exact intersection `[1.5, 3.5]`.

### Range divided by Function

[`divide_range_by_function()`](../../packages/treetime-distribution/src/distribution_ops/divide.rs) uses the same knot-filtering approach. A partial overlap can shrink to interior knots or become a `Point` merely because only one stored knot lies inside a nonzero-width analytical intersection.

### Function divided by Function

[`divide_function_by_function()`](../../packages/treetime-distribution/src/distribution_ops/divide.rs) always uses the dividend domain. When the divisor has narrower support, `resample_range_dx()` evaluates it beyond that support using `GridFn` constant extrapolation. The result therefore extends outside the exact support intersection.

`Function x Function` and `Formula x Function` multiplication already compute exact analytical endpoints before constructing a uniform result grid. They are not part of this boundary defect.

## Required behavior

For each affected operation, compute the strict intersection
`[max(a_min, b_min), min(a_max, b_max)]` before evaluating either operand. An endpoint-only contact remains empty, matching the current multiplication semantics. The result grid must include both exact endpoints and use spacing no coarser than the finest relevant input.

## Related ticket

- [kb/tickets/distribution-preserve-analytic-overlap-boundaries.md](../tickets/distribution-preserve-analytic-overlap-boundaries.md)

## Related issue

- [M-distribution-support-boundary-semantics-unresolved.md](M-distribution-support-boundary-semantics-unresolved.md) - exact finite-support intersections can be fixed independently of the generic extrapolation policy
