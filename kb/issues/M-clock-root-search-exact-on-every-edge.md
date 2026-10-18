# Clock root search optimizes the split only next to the best node, and outliers leak into the fit

Phase 0a of [timetree-coalescent-prior-and-clock-rate-streamlining.md](../proposals/timetree-coalescent-prior-and-clock-rate-streamlining.md).

## Problem

- `fn find_best_root()` ([`find_best_root.rs`](../../packages/treetime/src/clock/find_best_root/find_best_root.rs))
  scores every node's `ClockSet`, then optimizes the root position along a branch only on the parent and child edges of
  the best node, by an 11-point grid by default (Brent and golden section on request). A better position on any other
  edge is never seen.
- With the default variance model (zero on internal edges, unit weight per tip) the regression is OLS, and the score
  of a root at fraction $x$ along an edge is exactly a quadratic in $x$: shifting the two edge messages moves
  $\sum d$, $\sum td$ linearly and $\sum d^2$ quadratically, while $\sum t$, $\sum t^2$ and the norm stay. Both
  objectives (estimated rate: $\chi^2$; fixed rate: $\chi^2$ at that rate) are quadratic, and the regression rate is
  linear in $x$, so the positive-rate constraint is an interval.
- Outliers still enter the fit: an outlier tip contributes no sums to its own edge's `clock_to_parent`, but
  `clock_from_child` is built from its date in the backward pass and in the leaf-edge cost function
  ([H-clock-regression-includes-filtered-outliers.md](H-clock-regression-includes-filtered-outliers.md)).

## Fix direction

- Search every edge. Under the OLS variance model the optimum on an edge is the vertex of the parabola through the
  scores at $x = 0, \frac12, 1$, clamped to $[0, 1]$ and to the interval where the rate is positive when that is
  required. Other variance models (`--covariation`, the explicit variance flags) keep the configured optimizer, applied
  to every edge, until 0c retires them.
- An outlier tip contributes the empty `ClockSet` to its edge in both directions.

## Validation

- Analytic: the closed-form split matches a dense scan of the same cost; the parabola through three points of a
  quadratic is exact.
- The new search's $\chi^2$ is never above the old one's; results agree wherever the old search found the optimum.
- `test_clock_regression_outliers_are_left_out_of_the_fit` passes and is no longer ignored.
