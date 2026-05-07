# Clock command collapses date intervals to their mean

`assign_dates` at [packages/treetime/src/commands/clock/assign_dates.rs#L19-L23](../../packages/treetime/src/commands/clock/assign_dates.rs#L19-L23) collapses every `DateOrRange` to its scalar mean via `DateOrRange::mean` and stores only an `f64` in the node payload's `time` field. Date intervals (e.g., "2020.0-2020.5") lose their range information and are treated identically to point dates.

## Impact

The clock regression treats interval midpoints as exact observations, ignoring the uncertainty that the interval represents. For samples with wide date ranges (common in historical or environmental sequences), the regression underestimates the variance of the date contribution, producing overconfident clock rate estimates.

## v0 comparison

v0 `clock_filter_methods.py` uses `np.mean(node.raw_date_constraint)` in the same way, so this is a shared limitation rather than a v0/v1 divergence. The timetree pipeline in both v0 and v1 preserves the interval for `load_date_constraints`, but the clock pipeline itself does not use it.

## Affected code

- Date collapse: [packages/treetime/src/commands/clock/assign_dates.rs#L19-L23](../../packages/treetime/src/commands/clock/assign_dates.rs#L19-L23)
- Downstream consumer: [packages/treetime/src/commands/clock/clock_regression.rs#L54-L67](../../packages/treetime/src/commands/clock/clock_regression.rs#L54-L67)

## Fix

Incorporate interval width into the variance model: add a date-uncertainty term to `ClockParams::variance_offset_leaf` proportional to the interval width squared, so that wider date ranges contribute less weight to the regression.

## Related issues

- Source: [M-clock-date-interval-collapsed-to-mean.md](../issues/M-clock-date-interval-collapsed-to-mean.md) -- delete after full resolution
