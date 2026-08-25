# Date and range representation is inconsistent across v1

v1 encodes a calendar date as a floating-point year (`YYYY.F`, for example `2015.3`) and spans of such dates as interval structs, but the naming and typing have drifted into several competing conventions. The canonical target is fixed by [../decisions/datetime-numeric-date-naming.md](../decisions/datetime-numeric-date-naming.md) (numeric date, `DateNumeric` newtype, `DateRangeNumeric`); the code has not yet been consolidated onto it.

## Inconsistencies

1. **Three names for the numeric-date value.** The same `YYYY.F` value is called `numdate` in augur JSON I/O ([packages/treetime/src/commands/timetree/output/augur_node_data.rs](../../packages/treetime/src/commands/timetree/output/augur_node_data.rs)), `CalendarTime` in the coalescent ([packages/treetime/src/coalescent/time_coordinate.rs](../../packages/treetime/src/coalescent/time_coordinate.rs)), and `year_fraction` in the datetime utilities ([packages/treetime-utils/src/datetime/year_fraction.rs](../../packages/treetime-utils/src/datetime/year_fraction.rs)). `year_fraction` is inaccurate: the year fraction is only the fractional part `F` (the `0.3`), not the whole `2015.3`. `CalendarTime` is ambiguous, because "calendar time" reads as an actual calendar date rather than the decimal-year encoding.

2. **`DateRange` name collision.** The identifier names two different types: `DateRange { begin: DateTime<Utc>, end: DateTime<Utc> }`, real calendar timestamps ([packages/treetime-utils/src/datetime/date_range.rs](../../packages/treetime-utils/src/datetime/date_range.rs)); and `DateRange { start: f64, end: f64 }`, a numeric-date span ([packages/treetime-io/src/dates_csv.rs](../../packages/treetime-io/src/dates_csv.rs)). A bare `DateRange` cannot be understood without tracing the import, and the two cannot coexist in one scope without qualification.

3. **Range field naming split.** Interval endpoints use three vocabularies: `start` / `end` (dominant: the `f64` `DateRange`, `StartEnd` and annotation segments in [packages/treetime-io/src/auspice_types.rs](../../packages/treetime-io/src/auspice_types.rs), GFF segments in [packages/treetime-io/src/gff.rs](../../packages/treetime-io/src/gff.rs)), `begin` / `end` (timestamp `DateRange`), and `time_min` / `time_max` (private `TimeRange`). `start` / `end` dominates by roughly ten to one.

4. **Local numeric-date interval duplicates.** Two local structs reimplement a numeric-date interval that belongs in the shared datetime layer. `branch_length_likelihood.rs` defines a private `struct TimeRange { time_min: f64, time_max: f64 }` ([packages/treetime/src/timetree/inference/branch_length_likelihood.rs](../../packages/treetime/src/timetree/inference/branch_length_likelihood.rs)), a numeric-date interval with a third field convention. The coalescent output model defines `struct SegmentInterval { start: f64, end: f64 }` ([packages/treetime/src/commands/timetree/output/coalescent.rs](../../packages/treetime/src/commands/timetree/output/coalescent.rs)) for skyline segment spans; it follows the dominant `start` / `end` naming but is still a one-off numeric-date interval awaiting unification.

## Impact

One concept reads several ways, so call sites share no vocabulary, the same name denotes different representations, and reuse is hidden behind local one-off types. No demonstrated runtime effect: the defects are naming, typing, and maintainability.

## Proposed solution

Consolidate onto the canonical types per the decision:

- Add `DateNumeric` as a newtype over `f64` (`#[serde(transparent)]`, so it serializes as a bare number) and `DateRangeNumeric { start: DateNumeric, end: DateNumeric }` in the `treetime-utils` datetime module. `DateNumeric` provides the operations the current `CalendarTime` exposes (`value`, `max`, `is_finite`).
- Replace `year_fraction` (for the whole value) and `CalendarTime` with `DateNumeric`. Keep functions that operate on the fractional part `F` named for the fraction.
- Replace the `f64` `DateRange`, the private `TimeRange`, and the coalescent `SegmentInterval` with `DateRangeNumeric`, which also resolves the name collision (the timestamp `DateRange` keeps its name).
- Rename the timestamp `DateRange` fields `begin` / `end` to `start` / `end`.
- Keep `numdate` / `num_date` only as external field names in the augur/Auspice JSON schema, which is an external contract and is not renamed.

Nothing is shipped, so the renames and field changes are free.

## Related issues

- Decision: [../decisions/datetime-numeric-date-naming.md](../decisions/datetime-numeric-date-naming.md)
