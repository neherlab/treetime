# date_uncertainty_due_to_rate default interval typo

## v0 location

`ClockTree.date_uncertainty_due_to_rate()` (`#ClockTree`, `#date_uncertainty_due_to_rate`) [packages/legacy/treetime/treetime/clock_tree.py#L1068](../../packages/legacy/treetime/treetime/clock_tree.py#L1068)

## Erratum

Default parameter `interval=(0.05, 0.095)` is a typo for `(0.05, 0.95)`. The second value should represent the 95th percentile, not 9.5th.

## Evidence

- All internal callers pass explicit `interval` arguments with the correct `(0.05, 0.95)`:
  - `ClockTree.get_confidence_interval()` (`#get_confidence_interval`) at [packages/legacy/treetime/treetime/clock_tree.py#L1128](../../packages/legacy/treetime/treetime/clock_tree.py#L1128) forwards its own correctly defaulted `interval=(0.05, 0.95)`
  - `ClockTree.get_max_posterior_region()` (`#get_max_posterior_region`) at [packages/legacy/treetime/treetime/clock_tree.py#L1180](../../packages/legacy/treetime/treetime/clock_tree.py#L1180) passes a computed interval
- The method is public (no `_` prefix), so the buggy default is reachable by external consumers

## v0 impact

Dead code path within treetime's own codebase. External library consumers calling `date_uncertainty_due_to_rate(node)` without `interval` would get a near-zero confidence interval (0.05 to 0.095 quantiles).

## v1 status

v1 uses `(0.05, 0.95)` correctly.
