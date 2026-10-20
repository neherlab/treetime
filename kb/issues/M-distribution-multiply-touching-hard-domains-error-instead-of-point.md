# Multiplying two ranges whose hard domains touch at one point raises an internal error

`fn multiply_range_range()` decides an empty result from `overlap_start >= overlap_end`, so exact endpoint contact (`overlap_start == overlap_end`) takes the empty path. The empty-result guard treats touching hard domains as _not_ disjoint (matching the `endpoint_contact` case of `test_multiply_hard_domains_disjoint`), so the operation raises an internal error instead of returning a result.

## Symptom and reproduction

Multiply two ranges whose supports meet at a single point, for example `range(1.0, 2.0) * range(2.0, 3.0)`, such as two date ranges that share an endpoint. `multiply_range_range` computes `overlap_start == overlap_end == 2.0`, takes the `overlap_start >= overlap_end` branch, and routes to `guarded_empty_result`. The hard domains touch, so `hard_domains_disjoint` returns `false` and the multiplication returns an internal error.

## Impact and scope

- Affects `fn multiply_range_range()` [packages/treetime-distribution/src/distribution_ops/multiply.rs#L111-L126](../../packages/treetime-distribution/src/distribution_ops/multiply.rs#L111-L126)
- Triggers only on exact floating-point endpoint contact, so it is rare in practice
- Both `Empty` and the internal error are wrong for a measure-zero overlap

## Root cause

The touching case is a degenerate overlap, not a disjoint pair. All other multiplication arms of two extended operands already handle it: `fn multiplication_support_intersection()` [packages/treetime-distribution/src/distribution_ops/multiply.rs#L346-L363](../../packages/treetime-distribution/src/distribution_ops/multiply.rs#L346-L363) returns `SupportIntersection::Point(t)` at endpoint contact, which produces a point mass. `multiply_range_function`, `multiply_functions`, `multiply_formula_formula`, `multiply_formula_function`, and `multiply_formula_range` use it (see `test_multiply_range_function_endpoint_contact_returns_point` and `test_multiply_function_function_endpoint_contact_returns_point`). The range/range arm does not use `SupportIntersection`, so it never takes the point-mass path.

## Fix approach

Use `multiplication_support_intersection()` in `multiply_range_range()`, as the sibling arms do: at `overlap_start == overlap_end`, return a point mass at that coordinate (amplitude from the product of the two amplitudes) instead of routing to the empty guard. Reserve the empty guard for `overlap_start > overlap_end` (disjoint). Add an endpoint-contact test for two ranges.
