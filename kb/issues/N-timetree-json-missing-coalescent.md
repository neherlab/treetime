# timetree.json missing coalescent and skyline parameters

`timetree.json` always contains only `clock_rate`, `intercept`, and `stats`
regardless of whether coalescent or skyline models were used. No Tc value,
skyline grid, or coalescent likelihood is recorded.

## Root cause

`write_clock_model()` at
[`clock_output.rs`](../../packages/treetime/src/commands/clock/clock_output.rs)
serializes only the `ClockModel` fields. Coalescent Tc and skyline results are
local variables in `run_timetree_estimation()` and are never passed to the
output writer.

v0 records the optimized Tc in `molecular_clock.txt` and the coal_LH in
`trace_run.log`.

## Related issues

- [Missing skyline output files](N-timetree-missing-skyline-output.md)
  skyline results not written to any file
- [Missing output files compared to v0](N-timetree-missing-output-files.md)
  umbrella issue
