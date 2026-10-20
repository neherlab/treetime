# Date column header matching breaks on hash

`fn normalize_csv_headers()` strips leading and trailing `#` from TSV/CSV column headers ([packages/treetime-io/src/csv.rs#L242-L246](../../packages/treetime-io/src/csv.rs#L242-L246)). `fn get_col_name()` ([packages/treetime-io/src/csv.rs#L138-L166](../../packages/treetime-io/src/csv.rs#L138-L166)) then compares the requested column names with the normalized headers, but does not normalize the requested names. When the metadata file uses `#name` as the column header and the user requests `#name`, the stripped header (`name`) no longer matches. The zika_20 dataset uses this pattern ([data/zika/20/metadata.tsv](../../data/zika/20/metadata.tsv)).

Both branches of `fn get_col_name()` have the mismatch:

- **CLI**: `--name-column` (alias of `--metadata-id-columns`) gives a list of candidate names, which the candidate branch compares case-insensitively ([packages/treetime-io/src/csv.rs#L152-L165](../../packages/treetime-io/src/csv.rs#L152-L165)). `treetime timetree` on `data/zika/20` with `--name-column '#name'` fails:

  ```
  When reading dates
  When reading file 'data/zika/20/metadata.tsv'
  Unable to find column:
    Looking for: #name
    Available columns are: name, date, country
  ```

  v0 accepts `--name-column '#name'` on the same file.
- **Golden-master runner**: the test support passes the fixture value as a provided name ([packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_support.rs#L90-L100](../../packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_support.rs#L90-L100)), which the exact-match branch compares ([packages/treetime-io/src/csv.rs#L143-L150](../../packages/treetime-io/src/csv.rs#L143-L150)). `zika_20` in [gm_runner_inputs.json](../../packages/treetime/src/timetree/inference/__tests__/__fixtures__/gm_runner_inputs.json) sets `"name_column": "#name"`, so the `#[case::zika_20("zika_20")]` lines are commented out under `// TODO:` notes in all three runner test files (`test_gm_runner_marginal_dense.rs`, `test_gm_runner_marginal_sparse.rs`, `test_gm_runner_poisson.rs`). No test exercises the mismatch

Note: the `#`-stripping itself is a v1 improvement over v0. v0 does not strip `#` from column headers, so `--name-column=name` fails on zika/20 metadata (header is `#name`). v1 strips the `#` and finds "name" automatically. The bug is specifically the mismatch when the user passes `--name-column="#name"` (pre-stripped value does not match post-stripped header).

## Fix direction

Normalize the requested column names in the same way as the headers, then re-enable the commented-out `zika_20` cases of the runner tests.
