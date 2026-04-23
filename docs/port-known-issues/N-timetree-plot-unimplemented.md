# --plot-rtt and --plot-tree return error

Both `--plot-rtt` and `--plot-tree` flags are accepted by clap but return an
error at runtime with "not yet implemented" via `make_error!()`.

## Location

[`run.rs#L503-L509`](../../packages/treetime/src/commands/timetree/run.rs#L503-L509)

## Related issues

- N-timetree-plot-arg-type (resolved): the flags previously had the wrong type (`Option<usize>` instead of `Option<PathBuf>`)
- [Missing output files compared to v0](N-timetree-missing-output-files.md)
  `root_to_tip_regression.pdf` and `timetree.pdf` listed as missing outputs
