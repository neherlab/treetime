# Clock command collapses date intervals to their mean

> [!IMPORTANT]
> **Decision required.** v0 also regresses on the interval mean (see v0 comparison), so weighting intervals by their width would diverge from v0 and needs approval. Options: (a) keep v0 parity and treat the interval midpoint as an exact observation; (b) add a date-uncertainty term to the leaf variance, for example proportional to the squared interval width, so that wider date ranges carry less weight in the regression.

`assign_dates` at [packages/treetime/src/clock/assign_dates.rs#L28-L32](../../packages/treetime/src/clock/assign_dates.rs#L28-L32) collapses every `DateOrRange` to its scalar mean via `DateOrRange::mean` and stores only an `f64` in the node's `time` field. Date intervals (e.g., "2020.0-2020.5") lose their range information and are treated identically to point dates.

## Impact

The clock regression treats interval midpoints as exact observations, ignoring the uncertainty that the interval represents. For samples with wide date ranges (common in historical or environmental sequences), the regression underestimates the variance of the date contribution, producing overconfident clock rate estimates.

## v0 comparison

v0 uses `np.mean(node.raw_date_constraint)` in the same way, in the clock regression [packages/legacy/treetime/treetime/clock_tree.py#L275](../../packages/legacy/treetime/treetime/clock_tree.py#L275) and in the outlier filter [packages/legacy/treetime/treetime/clock_filter_methods.py#L12](../../packages/legacy/treetime/treetime/clock_filter_methods.py#L12), so this is a shared limitation rather than a v0/v1 divergence. The timetree pipeline in both v0 and v1 preserves the interval for `load_date_constraints`, but the clock pipeline itself does not use it.

A public support thread documents mixed-granularity sampling dates and recommends encoding unknown month or day components as `XX` [[issue](https://github.com/neherlab/treetime/issues/59)] [[comment](https://github.com/neherlab/treetime/issues/59#issuecomment-414949406)]. It establishes the user workflow that produces date intervals; it does not report the clock-regression weighting defect directly.

## Affected code

- Date collapse: [packages/treetime/src/clock/assign_dates.rs#L28-L32](../../packages/treetime/src/clock/assign_dates.rs#L28-L32)
- Downstream consumer: leaf date and leaf variance in the backward regression [packages/treetime/src/clock/clock_regression.rs#L341-L369](../../packages/treetime/src/clock/clock_regression.rs#L341-L369)
- Leaf variance parameter: `ClockVarianceParams::variance_offset_leaf` [packages/treetime/src/clock/clock_regression.rs#L452-L461](../../packages/treetime/src/clock/clock_regression.rs#L452-L461)

## Possible fix

If option (b) is approved: incorporate interval width into the variance model by adding a per-leaf date-uncertainty term to the leaf variance, proportional to the interval width squared, so that wider date ranges contribute less weight to the regression.
