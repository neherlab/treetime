# Timetree confidence interval is widened to contain the point estimate

> [!IMPORTANT]
> **Decision required.** v1 widens every confidence interval so that it contains the reported date. v0 and augur have no such step, so the step changes `num_date_confidence` against augur whenever the point estimate lies outside the computed region. Options:
>
> - Remove the step, for parity with v0 and augur
> - Keep the interval as computed and report the nodes whose date lies outside it, for example as a run warning
> - Keep the step and record it as an intentional change in `kb/decisions/`

## Summary

Confidence intervals of `timetree` combine the maximum posterior region of the marginal time posterior with the rate-susceptibility term of `--confidence`, in quadrature, clipped to the hard sides of the node's subtree distribution. They are written when `--confidence` is set.

## Details

`fn node_confidence_interval()` ([packages/treetime/src/timetree/confidence.rs#L114-L137](../../packages/treetime/src/timetree/confidence.rs#L114-L137)) takes the bounds from `fn confidence_bounds()` and then widens them with `lower.min(date)` and `upper.max(date)` ([packages/treetime/src/timetree/confidence.rs#L127-L128](../../packages/treetime/src/timetree/confidence.rs#L127-L128)). The step does not move the point estimate. A case where the computed interval and the reported date disagree (for example the clip to the hard sides of the node's subtree distribution) is hidden instead of reported.

v0 has no such step: `get_confidence_interval` and `get_max_posterior_region` return the `combine_confidence` result directly ([packages/legacy/treetime/treetime/clock_tree.py#L1142-L1144](../../packages/legacy/treetime/treetime/clock_tree.py#L1142-L1144), [packages/legacy/treetime/treetime/clock_tree.py#L1228-L1230](../../packages/legacy/treetime/treetime/clock_tree.py#L1228-L1230)).

### Effect on `num_date_confidence`

`num_date_confidence` in `timetree.node-data.json` (and `num_date.confidence` in `auspice_tree.json`) is the `[lower, upper]` 90% region of this module. Augur sets `num_date_confidence = list(tt.get_max_posterior_region(n, 0.9))`, and v1 composes the region the same way. A node with a point posterior (for example a tip with an exact date) gets `[date, date]`, as augur does through v0 `get_max_posterior_region` ([packages/legacy/treetime/treetime/clock_tree.py#L1170-L1171](../../packages/legacy/treetime/treetime/clock_tree.py#L1170-L1171)). The bounds differ from augur only where the widening step moves a bound to the point estimate, and on the posteriors of [kb/v0-errata/timetree-max-posterior-region-multimodal-side.md](../v0-errata/timetree-max-posterior-region-multimodal-side.md) and [kb/v0-errata/timetree-max-posterior-region-brent-plateau.md](../v0-errata/timetree-max-posterior-region-brent-plateau.md). Auspice colors and HPD bars read this field, so a change of the widening step changes the Auspice output as well as the node data file.

## Impact

- Inconsistencies between the interval and the point estimate are hidden instead of reported
- `num_date_confidence` differs from augur where the widening step applies
