# Column auto-detection gaps in CSV readers

[`get_col_name()`](../../packages/treetime-io/src/csv.rs#L125) (`#get_col_name`) is the shared column detection function used by both date and discrete state CSV readers. Three behavioral gaps exist relative to v0.

## D1. Date column: exact match vs substring

v1 [`read_dates_from_reader()`](../../packages/treetime-io/src/dates_csv.rs#L55) (`#read_dates_from_reader`) passes `["date"]` as the candidate list, requiring an exact header match. v0 [`parse_dates()`](../../packages/legacy/treetime/treetime/utils.py#L268) (`#parse_dates`) uses substring matching: any column whose lowered name contains `"date"` qualifies (e.g. `num_date`, `collection_date`, `date_submitted`).

A metadata file with only a `num_date` column (no `date` column) would auto-detect in v0 but fail in v1 unless the user passes `--date-column=num_date`.

### Affected commands

- `clock` - [`run()`](../../packages/treetime/src/commands/clock/run.rs#L56) passes `name_column` and `date_column` to `read_dates()`
- `timetree` - [`initialization.rs`](../../packages/treetime/src/commands/timetree/initialization.rs#L59) passes `name_column` and `date_column` to `read_dates()`

### Project datasets

All datasets in `data/` have an exact `date` column, so this gap does not affect project test data. The gap affects external metadata files that use alternative date column names.

### Fix direction

Change the date column auto-detection to use case-insensitive substring matching. This requires either:

- A matching mode parameter on `get_col_name()` (e.g. `MatchMode::Exact` vs `MatchMode::ContainsIgnoreCase`)
- A separate detection path in `read_dates_from_reader()` that bypasses `get_col_name()` for date columns

Priority order when multiple columns contain `"date"`: leftmost column wins (v0 behavior, already the v1 iteration order).

## D2. Name column: case-sensitive match -- RESOLVED

Fixed in CLI args unification branch: `get_col_name()` auto-detect branch now lowercases both headers and candidates before comparison. Candidate list sourced from `default_name_candidates()` in `treetime-io/src/csv.rs`, propagated through `MetadataIdArgs`.

## D3. Mugration positional fallback

v0's mugration command at [`wrappers.py:847`](../../packages/legacy/treetime/treetime/wrappers.py#L847) falls back to `columns[0]` for the name column when no named match is found. It falls back to `columns[1]` for the attribute column at [line 858](../../packages/legacy/treetime/treetime/wrappers.py#L858). v1 [`read_discrete_attrs()`](../../packages/treetime-io/src/discrete_states_csv.rs#L32) returns an error for the name column (no match in candidates) and requires explicit `--attribute` for the value column (passes empty candidates `&[]` at [line 33](../../packages/treetime-io/src/discrete_states_csv.rs#L33)).

### Affected command

- `mugration` only

### Analysis

v1's behavior (fail with clear error) is defensible. Silently picking `columns[0]` as the name column is fragile: if the first column is not a taxon identifier, mugration produces wrong results without warning. v0 prints a message ("Using column 'X' as taxon name") but does not ask for confirmation.

### Approaches

- Keep v1 behavior: require explicit column specification when auto-detection fails. Users get an actionable error message listing available columns.
- Add positional fallback with warning: fall back to `columns[0]`/`columns[1]` but emit a warning. Matches v0 convenience at the cost of silent misuse risk.
- Add positional fallback for value column only: the name column has well-known candidates (`name`, `strain`, `accession`); if none match, an error is appropriate. The attribute column has no well-known candidates, so positional fallback to `columns[1]` is more defensible (v0 behavior). Add a warning.

### Recommendation

Keep v1 behavior for the name column (error on no match). For the attribute column, consider adding a positional fallback to `columns[1]` with a warning, matching v0 convenience for simple two-column files. This is lowest priority of the three gaps.

## D4. Name column: candidate-list priority vs header-position priority

v1 `get_col_name()` iterates the candidate list `["strain", "name", "accession"]` and returns the first candidate found in the headers. v0 `parse_dates()` at [`utils.py#L270-L298`](../../packages/legacy/treetime/treetime/utils.py#L270) collects all matching candidates as `(column_index, column_name)` tuples and calls `sorted()`, which sorts by column index -- the leftmost matching column in the CSV header wins.

When a metadata file has both `accession` (column 0, matching tree leaf names) and `strain` (column 2, not matching), v0 picks `accession` (leftmost) and succeeds. v1 picks `strain` (first in candidate list) and the names don't match tree leaves, causing zero or near-zero date matches and downstream failures.

### Affected commands

All commands using `get_col_name()` with `default_name_candidates()`:

- `clock` -- via `MetadataIdArgs` in [`clock/args.rs`](../../packages/treetime/src/commands/clock/args.rs)
- `timetree` -- via `MetadataIdArgs` in [`timetree/args.rs`](../../packages/treetime/src/commands/timetree/args.rs)
- `mugration` -- via `MetadataIdArgs` in [`mugration/args.rs`](../../packages/treetime/src/commands/mugration/args.rs)

### Affected datasets

Datasets where tree leaves use accession identifiers and metadata has both `accession` and `strain` columns:

- `lassa/L/*` -- all sizes (20, 50, 200)
- `rsv/a/*` -- all sizes (20, 100)
- `mpox/clade-ii/*` -- size 20

### Fix direction

Change `get_col_name()` auto-detection to match v0 behavior: when no `provided_name` is given, find all candidates present in the headers and pick the one with the lowest column index (leftmost in header). The implementation at [`csv.rs#L144-L161`](../../packages/treetime-io/src/csv.rs#L144) currently uses `.find_map()` over candidates (candidate-list order). Replace with iteration over headers checking membership in candidates (header-position order).

## Related issues

- [Date column header matching breaks on hash](M-timetree-date-header-hash.md) - related `#`-stripping interaction with `--name-column` argument
