# Formula discretization errors silently swallowed

When `discretize_formula()` fails on a `Formula` variant, several methods on [`Distribution`](../../packages/treetime-distribution/src/distribution_core/distribution.rs#L193-L338) (`#Distribution`) return sentinel values instead of propagating the error:

- `max_value()` returns `0.0`
- `scale_by()` returns `Distribution::Empty`
- `to_neglog()` returns `Distribution::Empty`

These sentinel values are indistinguishable from legitimate states (`Empty` distribution, zero amplitude), so callers cannot detect discretization failure. The affected methods do not return `Result`, making error propagation impossible without API changes.

`Formula` distributions are created by coalescent code (`coalescent_contribution_node()`, `coalescent_contribution_merger()`, `skyline_coalescent_contribution()`). If a coalescent formula fails to discretize, downstream `ScaledDistribution` operations and coalescent Tc optimization silently treat the distribution as empty.

`to_neglog()` has zero production callers and is not a practical concern.
