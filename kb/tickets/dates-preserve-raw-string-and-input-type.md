# Preserve raw date string and input type in date parsing

Replace `DateOrRange` with `DateConstraint` carrying the original input string and a richer `DateValue` enum that distinguishes exact dates, uncertain dates, and explicit ranges.

## Changes

1. Define new types in `packages/treetime-io/src/dates_csv.rs`:

```rust
pub struct DateExact { pub value: f64 }
pub struct DateRange { pub start: f64, pub end: f64 }

pub enum DateValue {
  Exact(DateExact),
  Uncertain(DateRange),
  Range(DateRange),
}

pub struct DateConstraint {
  pub raw: String,
  pub value: DateValue,
}
```

2. Update `read_date()` to return `Option<DateConstraint>`, preserving the trimmed input string and tagging each parse path with the correct variant.

3. Update `DatesMap` to `BTreeMap<String, Option<DateConstraint>>`.

4. Update `convert_record()` and `DateRecord` accordingly.

5. Update consumers:
   - `packages/treetime/src/clock/date_constraints.rs:14-19` - `date_or_range_to_distribution()` matches on new enum
   - `packages/treetime/src/clock/assign_dates.rs` - uses `DatesMap`

6. Add accessor methods on `DateValue`: `mean()`, `is_exact()`. Add `DateRange::width()`, `DateRange::contains()`.

7. Add `year_fraction_to_datestring(f64) -> String` utility in `packages/treetime-utils/src/datetime/year_fraction.rs` using existing `year_fraction_to_date()` + `format("%Y-%m-%d")`.

8. Update tests in `packages/treetime-io/src/__tests__/test_dates_csv.rs`.

9. Remove `DateOrRange` enum (fully replaced).

## Naming collision

Existing `DateRange` in `packages/treetime-utils/src/datetime/date_range.rs` wraps `(DateTime<Utc>, DateTime<Utc>)`. The new `DateRange` wraps `(f64, f64)` year fractions. These live in different modules (`treetime-io` vs `treetime-utils`) and serve different purposes. If collision becomes an issue, the existing chrono-based type could be renamed to `CalendarDateRange`.

## Related issues

- Source: [M-dates-raw-string-not-preserved.md](../issues/M-dates-raw-string-not-preserved.md)
