# Clock command collapses date intervals to their mean

> [!IMPORTANT]
> **Decision required.** v0 also regresses on the interval mean (see v0 comparison), so weighting intervals by their width would diverge from v0 and needs approval. Options: (a) keep v0 parity and treat the interval midpoint as an exact observation; (b) weight each tip by its date uncertainty, for example inversely to the squared interval width, so that wider date ranges carry less weight in the regression. The regression is ordinary least squares with unit tip weights ([kb/decisions/clock-ordinary-least-squares-regression.md](../decisions/clock-ordinary-least-squares-regression.md)), so (b) would reintroduce per-tip weights in `ClockSet::leaf_contribution`.

`assign_dates` at [packages/treetime/src/clock/assign_dates.rs#L25-L28](../../packages/treetime/src/clock/assign_dates.rs#L25-L28) collapses every `DateConstraint` ([packages/treetime-primitives/src/date.rs](../../packages/treetime-primitives/src/date.rs)) to its scalar mean via `DateConstraint::mean` and stores only an `f64` in the node's `time` field. Date intervals (e.g., "2020.0-2020.5") lose their range information and are treated identically to point dates.

## Impact

The clock regression treats interval midpoints as exact observations with the same weight as exact dates, ignoring the uncertainty that the interval represents. For samples with wide date ranges (common in historical or environmental sequences), a midpoint far from the true sampling date pulls the estimated rate and root date as much as an exact date does. `clock` reports no rate uncertainty ([M-clock-command-reports-no-rate-uncertainty.md](M-clock-command-reports-no-rate-uncertainty.md)), so the effect is on the point estimates only.

## v0 comparison

v0 uses `np.mean(node.raw_date_constraint)` in the same way, in the clock regression [packages/legacy/treetime/treetime/clock_tree.py#L275](../../packages/legacy/treetime/treetime/clock_tree.py#L275) and in the outlier filter [packages/legacy/treetime/treetime/clock_filter_methods.py#L12](../../packages/legacy/treetime/treetime/clock_filter_methods.py#L12), so this is a shared limitation rather than a v0/v1 divergence. The timetree pipeline in both v0 and v1 preserves the interval for `load_date_constraints`, but the clock pipeline itself does not use it.

A public support thread documents mixed-granularity sampling dates and recommends encoding unknown month or day components as `XX` [[issue](https://github.com/neherlab/treetime/issues/59)] [[comment](https://github.com/neherlab/treetime/issues/59#issuecomment-414949406)]. It establishes the user workflow that produces date intervals; it does not report the clock-regression weighting defect directly.

## Affected code

- Date collapse: [packages/treetime/src/clock/assign_dates.rs#L25-L28](../../packages/treetime/src/clock/assign_dates.rs#L25-L28)
- Downstream consumer: leaf date in the backward regression, `fn ClockRegressionPass::backward_node()` [packages/treetime/src/clock/clock_regression.rs#L285](../../packages/treetime/src/clock/clock_regression.rs#L285), and the tip sums `fn ClockSet::leaf_contribution()` [packages/treetime/src/clock/clock_set.rs#L19](../../packages/treetime/src/clock/clock_set.rs#L19)

## Possible fix

If option (b) is approved: give each tip a weight that decreases with the squared interval width, carried in the tip sums of `ClockSet`, so that wider date ranges contribute less weight to the regression.
