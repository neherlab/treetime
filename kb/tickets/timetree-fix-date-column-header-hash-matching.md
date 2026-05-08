# Fix date column header matching when name contains hash

`read_dates()` (`#read_dates`) strips leading `#` from TSV/CSV column headers. When the metadata file uses `#name` as the column header and the CLI argument specifies `name_column="#name"`, the stripped header (`name`) no longer matches the expected column name (`#name`). The zika_20 dataset uses this pattern, and all three timetree runner test files have zika_20 commented out, including [test_gm_runner_poisson.rs#L26](../../packages/treetime/src/commands/timetree/inference/__tests__/test_gm_runner/test_gm_runner_poisson.rs#L26).

Note: the `#`-stripping itself is a v1 improvement over v0. v0 does not strip `#` from column headers, so `--name-column=name` fails on zika/20 metadata (header is `#name`). v1 strips the `#` and finds "name" automatically. The bug is specifically the mismatch when the user passes `--name-column="#name"` (pre-stripped value does not match post-stripped header).

## Related issues

- Source: [M-timetree-date-header-hash.md](../issues/M-timetree-date-header-hash.md) -- delete after full resolution
