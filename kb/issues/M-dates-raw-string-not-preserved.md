# Date parsing discards original input string

`read_date()` in `packages/treetime-io/src/dates_csv.rs:112-145` parses date strings into `DateOrRange` (numeric year fractions) and discards the original string. The `DateOrRange` enum also collapses distinct input types: `"2013-12-XX"` (uncertain) and `"2020-01-01/2020-06-30"` (explicit range) both become `YearFractionRange((f64, f64))`.

## Impact

Node data JSON output requires three fields derived from date input:

- `raw_date`: the original string from metadata (e.g. `"2013-12-XX"`)
- `date`: resolved date string in `"YYYY-MM-DD"` format (reverse of year fraction)
- `date_inferred`: whether the date was inferred vs known exactly

Without the original string, `raw_date` cannot be produced. Without distinguishing uncertain from exact inputs, `date_inferred` requires a separate tracking mechanism.

The reverse conversion (year fraction to `"YYYY-MM-DD"`) is available via `year_fraction_to_date()` in `packages/treetime-utils/src/datetime/year_fraction.rs:28-35`, so `date` is producible from `numdate`.

## Design

Replace `DateOrRange` with a richer type that preserves input semantics:

```rust
pub struct DateExact { pub value: f64 }
pub struct DateRange { pub start: f64, pub end: f64 }

pub enum DateValue {
  Exact(DateExact),       // "2013-12-31" or "2014.3"
  Uncertain(DateRange),  // "2013-12-XX" (ambiguous components)
  Range(DateRange),      // "2020-01-01/2020-06-30" (explicit range)
}

pub struct DateConstraint {
  pub raw: String,        // original input string
  pub value: DateValue,   // parsed numeric representation
}
```

`DateRange` is shared between `Uncertain` and `Range` to avoid duplicating range logic. `DateExact` is a newtype for future extensibility. `date_inferred` derives from the `DateValue` variant. `DatesMap` becomes `BTreeMap<String, Option<DateConstraint>>`.

## Locations

- `packages/treetime-io/src/dates_csv.rs:17-33` - `DateOrRange`, `DatesMap`, `DateRecord` type definitions
- `packages/treetime-io/src/dates_csv.rs:112-145` - `read_date()` where string is discarded
- `packages/treetime/src/clock/date_constraints.rs:14-19` - `date_or_range_to_distribution()` matches on `DateOrRange`
- `packages/treetime/src/clock/assign_dates.rs` - uses `DatesMap`

Blast radius: 3 production files + test files in `packages/treetime-io/src/__tests__/test_dates_csv.rs`.
