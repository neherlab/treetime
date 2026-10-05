# --plot-rtt and --plot-tree return error

Both `--plot-rtt` and `--plot-tree` flags of `timetree` are accepted by clap but return an error at argument conversion with "not yet implemented" via `make_error!()`.

A public enhancement request asks TreeTime to emit root-to-tip plots at multiple analysis points so date parsing, clock filtering, and failures can be diagnosed [[issue](https://github.com/neherlab/treetime/issues/228)]. The request is related output context; this issue specifically tracks accepted Rust CLI flags that fail at runtime.

## Location

[`packages/app-commands/src/commands/timetree/args.rs#L151-L156`](../../packages/app-commands/src/commands/timetree/args.rs#L151-L156)
