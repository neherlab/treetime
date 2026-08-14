# Multiplying operands whose hard domains touch at one point raises an internal error

Some multiplication arms decide an empty result from `overlap_start >= overlap_end`, so exact endpoint contact (`overlap_start == overlap_end`) takes the empty path. The empty-result guard treats touching hard domains as _not_ disjoint (matching `test_multiply_hard_domains_disjoint`'s `endpoint_contact` case), so the operation now raises an internal error instead of returning a result.

## Symptom and reproduction

Multiply two ranges (or a formula with a range or another formula) whose supports meet at a single point, for example `range(1.0, 2.0) * range(2.0, 3.0)`. `multiply_range_range` computes `overlap_start == overlap_end == 2.0`, takes the `overlap_start >= overlap_end` branch, and routes to `guarded_empty_result`. The hard domains touch, so `hard_domains_disjoint` returns `false` and the multiplication returns an internal error.

## Impact and scope

- Affects `multiply_range_range` [`packages/treetime-distribution/src/distribution_ops/multiply.rs#L119`](../../packages/treetime-distribution/src/distribution_ops/multiply.rs#L119), `multiply_formula_formula` [`multiply.rs#L194`](../../packages/treetime-distribution/src/distribution_ops/multiply.rs#L194), and `multiply_formula_range` [`multiply.rs#L270`](../../packages/treetime-distribution/src/distribution_ops/multiply.rs#L270).
- Triggers only on exact floating-point endpoint contact, so it is rare in practice.
- Before the empty-result guard these arms silently returned `Empty` at touching; both `Empty` and the internal error are wrong for a measure-zero overlap.

## Root cause

The touching case is a degenerate overlap, not a genuinely disjoint pair. The sibling arms `multiply_range_function` and `multiply_function_function` already handle it correctly: `distribution_support_intersection` returns `SupportIntersection::Point(t)` at endpoint contact, which produces a point mass (see `test_multiply_range_function_endpoint_contact_returns_point` and `test_multiply_function_function_endpoint_contact_returns_point`). The range/range and formula arms do not use `SupportIntersection`, so they never take the point-mass path.

## Fix approach

Give the range/range and formula arms the same endpoint-contact handling as the function arms: at `overlap_start == overlap_end`, return a point mass at that coordinate (amplitude from the pointwise product) instead of routing to the empty guard. Reserve the empty guard for `overlap_start > overlap_end` (genuinely disjoint). This keeps the empty invariant intact and matches the point-mass behavior of the sibling arms.
