# Date and range representation is inconsistent across v1

v1 encodes a calendar date as a floating-point year (`YYYY.F`, for example `2015.3`) and spans of such dates as interval structs, but the naming and typing have drifted into several competing conventions. The canonical term is numeric date, with a `DateNumeric` newtype and a `DateRangeNumeric` interval; the code has not yet been consolidated onto it.

Numeric date is the term of v0 (`numdate`, `numeric_date()`), augur, BEAST, and the Auspice schema (`num_date`). The R ecosystem says decimal date (`lubridate::decimal_date()`, `ape`), but TreeTime reads and writes the augur and Auspice formats, so a second term at that boundary would bring back the inconsistency.

## Inconsistencies

1. **Two names for the numeric-date value.** The same `YYYY.F` value is called `numdate` in augur JSON I/O ([packages/app-output/src/augur_node_data_refine.rs](../../packages/app-output/src/augur_node_data_refine.rs)) and `year_fraction` in the datetime utilities ([packages/treetime-utils/src/datetime/year_fraction.rs](../../packages/treetime-utils/src/datetime/year_fraction.rs)). `year_fraction` is inaccurate: the year fraction is only the fractional part `F` (the `0.3`), not the whole `2015.3`.

2. **`DateRange` name collision.** The identifier names two different types: `DateRange { begin: DateTime<Utc>, end: DateTime<Utc> }`, real calendar timestamps ([packages/treetime-utils/src/datetime/date_range.rs#L10-L13](../../packages/treetime-utils/src/datetime/date_range.rs#L10-L13)); and `DateRange { start: f64, end: f64 }`, a numeric-date span ([packages/treetime-primitives/src/date.rs#L58-L62](../../packages/treetime-primitives/src/date.rs#L58-L62)). The same file defines `DateExact { value: f64 }` ([packages/treetime-primitives/src/date.rs#L53-L56](../../packages/treetime-primitives/src/date.rs#L53-L56)), a further wrapper of a single numeric date. A bare `DateRange` cannot be understood without tracing the import, and the two cannot coexist in one scope without qualification.

3. **Range field naming split.** Interval endpoints use three vocabularies: `start` / `end` (dominant: the `f64` `DateRange`, `StartEnd` and annotation segments in [packages/treetime-io/src/auspice_types.rs](../../packages/treetime-io/src/auspice_types.rs), GFF segments in [packages/treetime-io/src/gff.rs](../../packages/treetime-io/src/gff.rs)), `begin` / `end` (timestamp `DateRange`), and `time_min` / `time_max` (private `TimeRange`). `start` / `end` dominates by roughly ten to one.

4. **Local numeric-date interval duplicates.** Two local structs reimplement a numeric-date interval that belongs in the shared datetime layer. `branch_length_likelihood.rs` defines a private `struct TimeRange { time_min: f64, time_max: f64 }` ([packages/treetime/src/timetree/inference/branch_length_likelihood.rs](../../packages/treetime/src/timetree/inference/branch_length_likelihood.rs)), a numeric-date interval with a third field convention. The coalescent output model defines `struct SegmentInterval { start: f64, end: f64 }` ([packages/treetime/src/timetree/coalescent.rs](../../packages/treetime/src/timetree/coalescent.rs)) for skyline segment spans; it follows the dominant `start` / `end` naming but is still a one-off numeric-date interval awaiting unification.

## Impact

One concept reads several ways, so call sites share no vocabulary, the same name denotes different representations, and reuse is hidden behind local one-off types. No demonstrated runtime effect: the defects are naming, typing, and maintainability.

## Proposed solution

Consolidate onto the canonical types:

- Add `DateNumeric` as a newtype over `f64` (a newtype struct, which deser serializes as a bare number) and `DateRangeNumeric { start: DateNumeric, end: DateNumeric }` in the `treetime-utils` datetime module. `DateNumeric` provides `value`, `max`, and `is_finite`.
- Replace `year_fraction` (for the whole value) with `DateNumeric`. Keep functions that operate on the fractional part `F` named for the fraction.
- Replace `DateExact` with `DateNumeric`.
- Replace the `f64` `DateRange`, the private `TimeRange`, and the coalescent `SegmentInterval` with `DateRangeNumeric`, which also resolves the name collision (the timestamp `DateRange` keeps its name).
- Rename the timestamp `DateRange` fields `begin` / `end` to `start` / `end`.
- Keep `numdate` / `num_date` only as external field names in the augur/Auspice JSON schema, which is an external contract and is not renamed.

Nothing is shipped, so the renames and field changes are free.
