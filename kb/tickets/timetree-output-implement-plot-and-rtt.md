# Implement --plot-rtt and --plot-tree

Both `--plot-rtt` and `--plot-tree` flags are accepted by clap but return an error at runtime with "not yet implemented" via `make_error!()`.

## Location

[`packages/treetime/src/commands/timetree/run.rs#L195-L201`](../../packages/treetime/src/commands/timetree/run.rs#L195-L201)

## Related issues

- Source: [N-timetree-plot-unimplemented.md](../issues/N-timetree-plot-unimplemented.md) -- delete after full resolution
- [kb/issues/N-timetree-missing-output-files.md](../issues/N-timetree-missing-output-files.md)
  `root_to_tip_regression.pdf` and `timetree.pdf` listed as missing outputs
