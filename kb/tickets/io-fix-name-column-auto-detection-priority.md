# Fix name-column auto-detection to use header-position priority

`get_col_name()` in [`packages/treetime-io/src/csv.rs#L131-L162`](../../packages/treetime-io/src/csv.rs#L131) picks the first match from the candidate list `["strain", "name", "accession"]` when auto-detecting the name column. v0 picks the candidate with the lowest column index in the CSV header (leftmost wins).

When a metadata file has `accession` at column 0 (matching tree leaf names) and `strain` at column 2 (not matching), v0 picks `accession` and succeeds. v1 picks `strain` and all downstream name matching fails -- zero date matches, negative clock rates, or "insufficient dated leaves" errors.

## What to change

In the auto-detection branch of `get_col_name()` (`csv.rs#L144-L161`), replace the `.find_map()` iteration over candidates with iteration over headers. For each header (in order), check whether it matches any candidate (case-insensitive). Return the first header that matches. This gives leftmost-in-header priority, matching v0.

Current code (candidate-list order):

```rust
headers_lower
  .iter()
  .find_map(|candidate| { ... })
```

Target behavior (header-position order):

```rust
headers_lower
  .iter()
  .enumerate()
  .find_map(|(idx, header)| {
    if candidates_lower.contains(header) { Some(idx) } else { None }
  })
```

## Smoke test config

Also add `--name-column=accession` to `dev/run-smoke-tests` for lassa, rsv, and mpox clock/timetree entries. This is belt-and-suspenders -- the code fix addresses the root cause, but explicit `--name-column` in the smoke tests documents the dataset convention and protects against regressions.

## Verification

After the fix, run smoke tests for the affected datasets without `--name-column`:

```bash
./dev/docker/run ./dev/dev r treetime -- clock --tree=data/lassa/L/20/tree.nwk --dates=data/lassa/L/20/metadata.tsv --outdir=tmp/verify/clock/lassa/L/20
./dev/docker/run ./dev/dev r treetime -- clock --tree=data/rsv/a/20/tree.nwk --dates=data/rsv/a/20/metadata.tsv --outdir=tmp/verify/clock/rsv/a/20
```

The lassa datasets should succeed (dates will match). RSV may still fail with negative clock rate -- that is a data characteristic (weak temporal signal), not a name-matching bug.

## Related issues

- Source: [kb/issues/M-dates-column-auto-detection-gaps.md](../issues/M-dates-column-auto-detection-gaps.md) (section D4) -- delete D4 section after resolution
