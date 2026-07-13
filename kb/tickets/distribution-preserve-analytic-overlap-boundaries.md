# Preserve analytical overlap boundaries in mixed-support distribution operations

`Range x Function`, `Range / Function`, and `Function / Function` currently lose exact support intersections. Range/function paths select existing function knots inside the range, while function/function division evaluates over the entire dividend domain and extrapolates a narrower divisor.

## Scope

Change only these operations:

- `multiply_range_function()` in [`packages/treetime-distribution/src/distribution_ops/multiply.rs`](../../packages/treetime-distribution/src/distribution_ops/multiply.rs);
- `divide_range_by_function()` in [`packages/treetime-distribution/src/distribution_ops/divide.rs`](../../packages/treetime-distribution/src/distribution_ops/divide.rs); and
- `divide_function_by_function()` in the same division module.

Do not change `Function x Function` multiplication, any `Formula` operation, generic `GridFn` extrapolation, `YAxisPolicy` semantics, or convolution. Do not introduce a grid-point cap.

## Implementation contract

Use one strict intersection helper for all three paths. Reuse or relocate the existing tuple-based `multiplication_eval_range()` logic so the helper returns
`[max(a_min, b_min), min(a_max, b_max)]` only when `overlap_min < overlap_max`. Endpoint-only contact remains `Empty`, matching current multiplication behavior.

For each non-empty intersection:

1. Choose the finest relevant input spacing:
   - function spacing for `Range x Function` and `Range / Function`;
   - `min(dividend.dx(), divisor.dx())` for `Function / Function`.
2. Compute `n_points = ceil((overlap_max - overlap_min) / dx) + 1`, with at least two points.
3. Build a uniform grid whose first and last coordinates are the exact analytical intersection endpoints. Recomputing the realized spacing from the endpoint range is expected.
4. Evaluate both operands on that grid and apply the existing `YAxisPolicy::multiply`, `YAxisPolicy::divide`, and `YAxisPolicy::safe_divisor` behavior.

Because the intersection width cannot exceed the domain width of the input supplying the finest spacing, this construction is bounded by that input's length over the intersection. No independent cap or approximation policy is needed.

Remove the epsilon-expanded knot filters and the single-retained-knot conversion to `Point`; a positive-width analytical overlap always produces a function with both endpoints.

## Tests

Add parameterized unit tests for each operation. Cases must distinguish:

- full containment;
- left and right partial overlap with boundaries between input knots;
- a positive overlap containing no pre-existing interior knot;
- disjoint supports; and
- endpoint-only contact.

For multiplication, cover both operand orders. Compare the whole result distribution, including exact bounds and sampled values. For division, verify that no output point lies outside the divisor support.

Add property tests for support intersection, operand-order symmetry of range/function multiplication, and multiplication/division round trips where the divisor is positive and the mathematical inverse is valid. Use project property assertion helpers for floating-point comparisons.

Add a v0 golden-master comparison for a non-grid-aligned partial overlap. The capture must call v0 `Distribution.multiply()` and `Distribution.divide()` from [`packages/legacy/treetime/treetime/distribution.py`](../../packages/legacy/treetime/treetime/distribution.py), representing the range as a two-endpoint constant distribution. Preserve full precision and cite that v0 source next to every captured expected value or fixture assertion. Do not generate expectations through v1.

Validate new tests by temporarily perturbing an expected endpoint and confirming a clear failure before restoring it.

## Acceptance criteria

- All three scoped operations return the exact strict support intersection.
- Result spacing is no coarser than the finest relevant input spacing.
- Positive-width overlaps do not collapse to `Point` because of knot placement.
- Existing function/function multiplication and Formula behavior remain unchanged.
- Unit, property, and v0 golden-master tests cover non-grid-aligned boundaries.
- The full formatter, linter/build, and test suite pass:

  ```bash
  ./dev/docker/run ./dev/dev f
  ./dev/docker/run ./dev/dev l
  ./dev/docker/run ./dev/dev t
  ```

## Related issues

- Source: [kb/issues/M-distribution-mixed-support-operations-lose-exact-intersection-boundaries.md](../issues/M-distribution-mixed-support-operations-lose-exact-intersection-boundaries.md)
