# Mugration column detection has no positional fallback

v0's mugration falls back to header position when named column detection fails: [`columns[0]`](../../packages/legacy/treetime/treetime/wrappers.py#L847) for the taxon-name column and [`columns[1]`](../../packages/legacy/treetime/treetime/wrappers.py#L858) for the attribute column. v1 [`discrete_attrs_read()`](../../packages/treetime-io/src/discrete_states_csv.rs#L30-L48) hard-errors instead:

```rust
let name_column_idx = get_col_name(&headers, name_candidates, name_column)?;
let value_column_idx = get_col_name(&headers, &[], value_column)?;
```

[`get_col_name()`](../../packages/treetime-io/src/csv.rs#L138) returns an error when no candidate matches. The value column is queried with an empty candidate list, so it errors unless a column name is given.

The CLI never reaches that error for the value column: `--attribute` is a required argument, and argument validation rejects a run without it before the metadata file is read ([packages/app-commands/src/commands/mugration/args.rs#L37-L60](../../packages/app-commands/src/commands/mugration/args.rs#L37-L60)). A positional fallback for the attribute column therefore also means making `--attribute` optional.

## Affected command

- `mugration` only

## Analysis

v1's behavior (fail with a clear error listing available columns) is defensible. Silently picking `columns[0]` as the name column is fragile: if the first column is not a taxon identifier, mugration produces wrong results without warning. v0 prints a message ("Using column 'X' as taxon name") but does not ask for confirmation.

## Approaches

- Keep v1 behavior: require explicit column specification when auto-detection fails. Users get an actionable error listing available columns.
- Positional fallback with warning: fall back to `columns[0]`/`columns[1]` but emit a warning. Matches v0 convenience at the cost of silent-misuse risk.
- Positional fallback for the value column only: the name column has well-known candidates (`name`, `strain`, `accession`), so an error on no match is appropriate. The attribute column has no well-known candidates, so a positional fallback to `columns[1]` with a warning is more defensible (v0 behavior).

## Decision

> [!IMPORTANT]
> **Decision required.** v1 requires an explicit column where v0 falls back to a position. The options are the three approaches above. One candidate combination keeps the error for the name column, which has well-known candidates, and adds a positional fallback with a warning for the attribute column, which matches v0 for simple two-column files and requires an optional `--attribute`.

## Related issues

- [N-io-sequence-name-matching-unreliable.md](N-io-sequence-name-matching-unreliable.md) - related name matching defects in metadata readers
