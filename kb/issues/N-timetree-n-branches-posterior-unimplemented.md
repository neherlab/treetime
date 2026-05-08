# --n-branches-posterior returns error

The `--n-branches-posterior` flag is accepted by clap but returns an error at runtime with "not yet implemented" via `make_error!()`.

## Location

[`run.rs#L177`](../../packages/treetime/src/commands/timetree/run.rs#L177)
