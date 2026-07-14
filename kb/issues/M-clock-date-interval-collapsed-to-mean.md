# Clock command collapses date intervals to their mean

`assign_dates` at [packages/treetime/src/clock/assign_dates.rs#L19-L23](../../packages/treetime/src/clock/assign_dates.rs#L19-L23) collapses every `DateOrRange` to its scalar mean via `DateOrRange::mean` and stores only an `f64` in the node payload's `time` field. Date intervals (e.g., "2020.0-2020.5") lose their range information and are treated identically to point dates.

## Impact

The clock regression treats interval midpoints as exact observations, ignoring the uncertainty that the interval represents. For samples with wide date ranges (common in historical or environmental sequences), the regression underestimates the variance of the date contribution, producing overconfident clock rate estimates.

## v0 comparison

v0 `clock_filter_methods.py` uses `np.mean(node.raw_date_constraint)` in the same way, so this is a shared limitation rather than a v0/v1 divergence. The timetree pipeline in both v0 and v1 preserves the interval for `load_date_constraints`, but the clock pipeline itself does not use it.

A public support thread documents mixed-granularity sampling dates and recommends encoding unknown month or day components as `XX` [[issue](https://github.com/neherlab/treetime/issues/59)] [[comment](https://github.com/neherlab/treetime/issues/59#issuecomment-414949406)]. It establishes the user workflow that produces date intervals; it does not report the clock-regression weighting defect directly.

## Affected code

- Date collapse: [packages/treetime/src/clock/assign_dates.rs#L19-L23](../../packages/treetime/src/clock/assign_dates.rs#L19-L23)
- Downstream consumer: [packages/treetime/src/clock/clock_regression.rs#L54-L67](../../packages/treetime/src/clock/clock_regression.rs#L54-L67)

## Fix

Incorporate interval width into the variance model: add a date-uncertainty term to `ClockParams::variance_offset_leaf` proportional to the interval width squared, so that wider date ranges contribute less weight to the regression.

## Related tickets

- [kb/tickets/clock-date-interval-collapsed-to-mean.md](../tickets/clock-date-interval-collapsed-to-mean.md)
