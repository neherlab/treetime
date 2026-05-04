# --plot-rtt and --plot-tree return error

Both `--plot-rtt` and `--plot-tree` flags are accepted by clap but return an
error at runtime with "not yet implemented" via `make_error!()`.

## Location

[`run.rs#L503-L508`](../../packages/treetime/src/commands/timetree/run.rs#L503-L508)

## Related issues

- [Missing output files compared to v0](N-timetree-missing-output-files.md)
  `root_to_tip_regression.pdf` and `timetree.pdf` listed as missing outputs
