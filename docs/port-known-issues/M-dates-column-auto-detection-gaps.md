# Column auto-detection gaps in CSV readers

[`get_col_name()`](../../packages/treetime-io/src/csv.rs#L125) (`#get_col_name`) is the shared column detection function used by both date and discrete state CSV readers. Three behavioral gaps exist relative to v0.

## D1. Date column: exact match vs substring

v1 [`read_dates_from_reader()`](../../packages/treetime-io/src/dates_csv.rs#L55) (`#read_dates_from_reader`) passes `["date"]` as the candidate list, requiring an exact header match. v0 [`parse_dates()`](../../packages/legacy/treetime/treetime/utils.py#L268) (`#parse_dates`) uses substring matching: any column whose lowered name contains `"date"` qualifies (e.g. `num_date`, `collection_date`, `date_submitted`).

A metadata file with only a `num_date` column (no `date` column) would auto-detect in v0 but fail in v1 unless the user passes `--date-column=num_date`.

### Affected commands

- `clock` - [`run()`](../../packages/treetime/src/commands/clock/run.rs#L64) passes `name_column` and `date_column` to `read_dates()`
- `timetree` - [`initialization.rs`](../../packages/treetime/src/commands/timetree/initialization.rs#L59) passes `name_column` and `date_column` to `read_dates()`

### Project datasets

All datasets in `data/` have an exact `date` column, so this gap does not affect project test data. The gap affects external metadata files that use alternative date column names.

### Fix direction

Change the date column auto-detection to use case-insensitive substring matching. This requires either:

- A matching mode parameter on `get_col_name()` (e.g. `MatchMode::Exact` vs `MatchMode::ContainsIgnoreCase`)
- A separate detection path in `read_dates_from_reader()` that bypasses `get_col_name()` for date columns

Priority order when multiple columns contain `"date"`: leftmost column wins (v0 behavior, already the v1 iteration order).

## D2. Name column: case-sensitive match

v1 `get_col_name()` at [line 141](../../packages/treetime-io/src/csv.rs#L141) uses `possible_names.contains(header)`, which is case-sensitive. v0 `parse_dates()` at [line 270](../../packages/legacy/treetime/treetime/utils.py#L270) lowercases the header before comparing: `col.lower()`. A header like `Name`, `STRAIN`, or `Accession` would auto-detect in v0 but fail in v1.

### Affected commands

- `clock` and `timetree` (via `read_dates()`)
- `mugration` (via [`read_discrete_attrs()`](../../packages/treetime-io/src/discrete_states_csv.rs#L32))

### Fix direction

Lowercase both the header and the candidate names before comparison in `get_col_name()`. The headers are already preprocessed (trimmed, `#`-stripped) in both `dates_csv.rs` and `discrete_states_csv.rs`, so adding `.to_lowercase()` to the comparison is the natural extension.

Note: v0's mugration path at [`wrappers.py:840-845`](../../packages/legacy/treetime/treetime/wrappers.py#L840) is case-sensitive (checks `"name"`, `"strain"`, `"accession"` against raw column names). v1 case-insensitive matching would be a superset of both v0 behaviors, which is acceptable.

## D3. Mugration positional fallback

v0's mugration command at [`wrappers.py:847`](../../packages/legacy/treetime/treetime/wrappers.py#L847) falls back to `columns[0]` for the name column when no named match is found. It falls back to `columns[1]` for the attribute column at [line 858](../../packages/legacy/treetime/treetime/wrappers.py#L858). v1 [`read_discrete_attrs()`](../../packages/treetime-io/src/discrete_states_csv.rs#L32) returns an error for the name column (no match in candidates) and requires explicit `--attribute` for the value column (passes empty candidates `&[]` at [line 33](../../packages/treetime-io/src/discrete_states_csv.rs#L33)).

### Affected command

- `mugration` only

### Analysis

v1's behavior (fail with clear error) is defensible. Silently picking `columns[0]` as the name column is fragile: if the first column is not a taxon identifier, mugration produces wrong results without warning. v0 prints a message ("Using column 'X' as taxon name") but does not ask for confirmation.

### Approaches

- **Keep v1 behavior**: require explicit column specification when auto-detection fails. Users get an actionable error message listing available columns.
- **Add positional fallback with warning**: fall back to `columns[0]`/`columns[1]` but emit a warning. Matches v0 convenience at the cost of silent misuse risk.
- **Add positional fallback for value column only**: the name column has well-known candidates (`name`, `strain`, `accession`); if none match, an error is appropriate. The attribute column has no well-known candidates, so positional fallback to `columns[1]` is more defensible (v0 behavior). Add a warning.

### Recommendation

Keep v1 behavior for the name column (error on no match). For the attribute column, consider adding a positional fallback to `columns[1]` with a warning, matching v0 convenience for simple two-column files. This is lowest priority of the three gaps.

## Related issues

- [Date column header matching breaks on hash](M-timetree-date-header-hash.md) - related `#`-stripping interaction with `--name-column` argument
