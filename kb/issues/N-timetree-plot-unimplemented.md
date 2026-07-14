# --plot-rtt and --plot-tree return error

Both `--plot-rtt` and `--plot-tree` flags are accepted by clap but return an error at runtime with "not yet implemented" via `make_error!()`.

A public enhancement request asks TreeTime to emit root-to-tip plots at multiple analysis points so date parsing, clock filtering, and failures can be diagnosed [[issue](https://github.com/neherlab/treetime/issues/228)]. The request is related output context; this issue specifically tracks accepted Rust CLI flags that fail at runtime.

## Location

[`packages/treetime/src/commands/timetree/run.rs#L195-L201`](../../packages/treetime/src/commands/timetree/run.rs#L195-L201)

## Related issues

- [Missing output files compared to v0](N-timetree-missing-output-files.md)
  `root_to_tip_regression.pdf` and `timetree.pdf` listed as missing outputs

## Related tickets

- [kb/tickets/timetree-output-implement-plot-and-rtt.md](../tickets/timetree-output-implement-plot-and-rtt.md)
