# --plot-rtt and --plot-tree panic with todo!()

Both `--plot-rtt` and `--plot-tree` flags are accepted by clap but crash at
runtime with "not yet implemented: plot_root_to_tip/plot_time_tree not yet
implemented".

## Location

[`run.rs#L258-L262`](../../packages/treetime/src/commands/timetree/run.rs#L258-L262)

## Related issues

- [--plot-\* argument typed as Option\<usize\>](N-timetree-plot-arg-type.md)
  the flags also have the wrong type (should be `Option<PathBuf>`)
- [Missing output files compared to v0](N-timetree-missing-output-files.md)
  `root_to_tip_regression.pdf` and `timetree.pdf` listed as missing outputs
