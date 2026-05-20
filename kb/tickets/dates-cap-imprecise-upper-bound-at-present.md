# Cap imprecise date upper bound at present

v0 [`ambiguous_date_to_date_range()`](../../packages/legacy/treetime/treetime/utils.py#L345) (`#ambiguous_date_to_date_range`) caps the upper bound of an uncertain date range at today's date (line 396). v1 [`determine_date_range()`](../../packages/treetime-utils/src/datetime/parse_uncertain_date.rs#L34) (`#determine_date_range`) returns the full calendar range without capping.

A date like `2026-XX-XX` produces range `[2026-01-01, 2026-12-31]` in v1. If parsed before December 31, v0 would cap the upper bound at today, producing a narrower range. v1 keeps the full year.

## Affected commands

- `clock` - [`assign_dates()`](../../packages/treetime/src/clock/assign_dates.rs#L9) (`#assign_dates`) uses `DateOrRange::mean()` for regression. A future upper bound shifts the midpoint forward.
- `timetree` - [`load_date_constraints()`](../../packages/treetime/src/clock/date_constraints.rs#L23) (`#load_date_constraints`) converts range to `Distribution::range()`. A future upper bound widens the uniform prior, allowing the optimizer to place the node in the future.

## Science background

Molecular clock inference assumes leaf sample dates are known. Imprecise dates (e.g. `2020-XX-XX`) represent collection date uncertainty, not future sampling. A sample collected in 2026 cannot have been collected on 2026-12-31 if today is 2026-03-15. Capping at today is scientifically correct for leaf dates.

For internal node date constraints (rare in practice), the argument is weaker. An internal node date constraint represents a calibration point, and the full calendar range may be appropriate.

## Approaches

### A1. Cap at today (v0 parity)

Cap the upper bound at today in [`read_date()`](../../packages/treetime-io/src/dates_csv.rs#L112) (`#read_date`) after parsing an uncertain date. Straightforward v0 parity. Introduces a non-determinism source (output depends on execution date).

### A2. Cap at today with warning

Same as A1, but emit a warning when capping occurs. Makes the behavior visible to users. Preferred over silent capping.

### A3. Hard error on future upper bound

Reject uncertain dates whose resolved upper bound exceeds today. Strictest option. Breaks workflows where users intentionally specify wide ranges spanning the present.

### A4. No cap, accept the difference

Keep v1 behavior. The full calendar range is the literal interpretation of the input. For dates far in the past (the common case), capping has no effect. For recent dates, the difference between `2026-03-15` and `2026-12-31` as an upper bound is small relative to typical date uncertainty in phylogenetic datasets.

### A5. Configurable via CLI flag

Add a flag (e.g. `--cap-uncertain-dates`) controlling behavior. Adds complexity for a rare edge case.

## Recommendation

A2 (cap with warning) matches v0 behavior, is scientifically defensible for leaf dates, and the warning makes the non-obvious behavior visible. The implementation point is [`read_date()`](../../packages/treetime-io/src/dates_csv.rs#L112) in `dates_csv.rs`, where `DateOrRange::YearFractionRange` values from uncertain parsing can be clamped before returning.

## Related issues

- Source: [N-dates-imprecise-upper-bound-not-capped.md](../issues/N-dates-imprecise-upper-bound-not-capped.md) -- delete after full resolution
