# Date column header matching breaks on hash

`read_dates()` (`#read_dates`) strips leading `#` from TSV/CSV column headers.
When the metadata file uses `#name` as the column header and the CLI argument
specifies `name_column="#name"`, the stripped header (`name`) no longer matches
the expected column name (`#name`). The zika_20 dataset uses this pattern, and
all three timetree runner test files have zika_20 commented out, including
[test_gm_runner_poisson.rs#L26](../../packages/treetime/src/commands/timetree/inference/__tests__/test_gm_runner/test_gm_runner_poisson.rs#L26).
