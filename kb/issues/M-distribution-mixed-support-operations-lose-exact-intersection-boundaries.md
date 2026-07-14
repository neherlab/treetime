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

## Endpoint-contact semantics are undecided

The current implementation treats exact endpoint contact as `Empty` because `multiplication_eval_range()` requires `eval_min < eval_max`. [`packages/treetime-distribution/src/distribution_ops/multiply.rs#L143-L162`](../../packages/treetime-distribution/src/distribution_ops/multiply.rs#L143-L162)

Let $\varepsilon$ denote v0's numerical tolerance. V0 consolidates operand knots within an $\varepsilon$-expanded intersection and converts a single surviving knot into a delta distribution. [`packages/legacy/treetime/treetime/distribution.py#L82-L145`](../../packages/legacy/treetime/treetime/distribution.py#L82-L145) Exact endpoint contact therefore produces a point distribution in v0 rather than `Empty`.

Discarding point contact is an unapproved parity divergence. Choose one contract before implementation:

- match v0 by returning a point distribution evaluated at the shared endpoint; or
- explicitly approve positive-measure support semantics, document why point contact is discarded, and define how explicit point distributions interact with that rule.

## Required behavior for positive-width intersections

Let $[a_{\min},a_{\max}]$ and $[b_{\min},b_{\max}]$ denote the operand supports. For each affected operation, compute the strict intersection $\left[\max(a_{\min},b_{\min}),\min(a_{\max},b_{\max})\right]$ before evaluating either operand. The result grid must include both exact endpoints and use spacing no coarser than the finest relevant input.

No implementation ticket is ready until endpoint-contact semantics are approved.

## Related issue

- [M-distribution-support-boundary-semantics-unresolved.md](M-distribution-support-boundary-semantics-unresolved.md) - exact finite-support intersections can be fixed independently of the generic extrapolation policy
