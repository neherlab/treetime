# --n-branches-posterior panics with todo!()

The `--n-branches-posterior` flag is accepted by clap but crashes at runtime with
"not yet implemented: n_branches_posterior".

## Location

[`run.rs#L113`](../../packages/treetime/src/commands/timetree/run.rs#L113)
