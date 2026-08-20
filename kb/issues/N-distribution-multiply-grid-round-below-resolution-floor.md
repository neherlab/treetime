# Multiplication/division grid rounding can fall below the finest-operand resolution floor

The tail-aware arithmetic grid contract sizes the result grid as
`n = clamp(round((x_max - x_min) / dx) + 1, 2, 1e6)` with `dx = min` of the operand spacings
(`distribution_support_n_points` in `packages/treetime-distribution/src/distribution_ops/time_bounds.rs`).
Because the point count is rounded to the nearest integer, the realized spacing `(x_max - x_min)/(n-1)`
can come out slightly coarser than the finest operand's `dx`.

## Symptom and reproduction

For an intersection width that is a non-integer multiple of `dx`, `round` picks the lower point count
when the fractional part is below `0.5`. Example: width `5.4 * dx` rounds to `5` intervals, giving a
realized spacing of `1.08 * dx`, coarser than the finest operand fed in.

## Impact and scope

- The proposal mandates a resolution floor: "The resolution floor is the finest operand `dx`"
  ([kb/proposals/distribution-log-space-and-hard-soft-boundaries-2.md](../proposals/distribution-log-space-and-hard-soft-boundaries-2.md), K23; also the "Point budget"
  in [kb/proposals/distribution-log-space-and-hard-soft-boundaries.md](../proposals/distribution-log-space-and-hard-soft-boundaries.md)). `round` can breach it.
- Applies to every arm that sizes a grid this way: `Function * Function`, `Range * Function`,
  `Formula * Function` (`multiply.rs`), and both division arms (`divide.rs`).
- Effect is sub-`dx` and bounded; no incorrect result, only a resolution slightly below target on some
  intersections. The endpoint-inclusion property is unaffected: the grid still lands both analytical
  endpoints on nodes via `Array1::linspace`, so a mode sitting on a hard bound stays representable.

## Current status

The decision [kb/decisions/distribution-tails-and-arithmetic.md](../decisions/distribution-tails-and-arithmetic.md) records the `round` formula and
explicitly lists this under "Accepted limitations" ("Rounding can make the realized spacing slightly
finer or coarser than the target"). So the landed behaviour is a documented, accepted limitation that
is now in tension with the later proposal's resolution-floor rule.

## Fix approach

Replace `round` with `ceil` in the point-count rule, so the realized spacing is always no coarser than
`dx` while both endpoints stay on grid nodes. This is a scientific-change-gated decision:

- It changes numerical output on non-integer-width intersections and must be applied consistently to
  multiplication and division, so the single grid formula in the decision does not fracture.
- It requires updating the consent-gated grid-contract formula in
  [kb/decisions/distribution-tails-and-arithmetic.md](../decisions/distribution-tails-and-arithmetic.md) and removing the matching "Accepted
  limitations" clause.
- The `round`-contract oracle tests (`test_multiply_function_function_uses_finer_spacing_over_intersection`,
  `test_multiply_formula_function_uses_function_spacing_over_intersection`) encode the current count and
  would need recomputed expectations under `ceil`.

Requires explicit human decision before implementation: adopt `ceil` per the proposal's resolution
floor, or keep `round` as the accepted limitation.
