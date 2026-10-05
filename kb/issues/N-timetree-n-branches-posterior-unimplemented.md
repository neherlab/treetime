# --n-branches-posterior returns error

The `--n-branches-posterior` flag is accepted by clap but returns an error at runtime with "not yet implemented" (`OperationError::InvalidParams` from the timetree pipeline).

## Location

[`pipeline.rs#L221-L223`](../../packages/treetime/src/timetree/pipeline.rs#L221-L223)
