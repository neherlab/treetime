# Date column auto-detection requires exact header match

When no `--date-column` is given, v1 auto-detects the date column by exact header equality, whereas v0 uses case-insensitive substring matching. A metadata file whose date column is named `num_date`, `collection_date`, or `date_submitted` auto-detects in v0 but fails in v1.

## Cause

[`read_metadata_table_from_reader()`](../../packages/treetime-io/src/dates_csv.rs#L114) passes the single candidate `["date"]` to the shared column detector [`get_col_name()`](../../packages/treetime-io/src/csv.rs#L55):

```rust
let date_index = get_col_name(&columns, &date_candidates, date_column);
```

The auto-detect branch of `get_col_name()` matches by whole-string equality after lowercasing, not substring containment:

```rust
candidates_lower.contains(&header_lower)
```

`Vec::contains` compares each candidate to the full header string, so `["date"]` matches a header named `date` (any case) but never `num_date`. v0 [`parse_dates()`](../../packages/legacy/treetime/treetime/utils.py#L268) instead qualifies any column whose lowered name _contains_ `"date"`.

## Affected commands

- `clock` - [`run.rs`](../../packages/app-commands/src/commands/clock/run.rs#L46) passes the ID column candidates and `date_column` to `read_dates()`
- `timetree` - [`initialization.rs`](../../packages/app-commands/src/commands/timetree/initialization.rs#L63) passes the ID column candidates and `date_column` to `read_dates()`
- `check-inputs` - [`check_inputs.rs`](../../packages/app-commands/src/check_inputs.rs) reads the same table through `read_metadata_table()`, so it reports no date column where the commands fail

The help text of `struct DateColumnArgs` in [`packages/app-commands/src/commands/shared/metadata.rs`](../../packages/app-commands/src/commands/shared/metadata.rs) describes the v0 rule ("the leftmost column whose name contains `date`"), so it disagrees with the current behavior.

## Project datasets

All datasets in `data/` have an exact `date` column, so this does not affect project test data. It affects external metadata files that use alternative date column names.

## Fix direction

Change date-column auto-detection to case-insensitive substring matching. Either:

- add a matching mode to `get_col_name()` (e.g. `MatchMode::Exact` vs `MatchMode::ContainsIgnoreCase`), or
- add a separate detection path in `read_metadata_table_from_reader()` that bypasses `get_col_name()` for date columns.

When multiple columns contain `"date"`, the leftmost column wins (v0 behavior, already the v1 header-iteration order in `get_col_name()`).

## Related issues

- [M-timetree-date-header-hash.md](M-timetree-date-header-hash.md) - `#`-stripping interaction with the `--name-column` argument
- [M-io-sequence-name-matching-unreliable.md](M-io-sequence-name-matching-unreliable.md) - related name matching defects in metadata readers
