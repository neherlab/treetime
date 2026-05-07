# Add log_lh parameter to DenseSeqDis::new()

## Description

`packages/treetime/src/representation/payload/dense.rs:53-55:`

`struct DenseSeqDis` holds a 2D probability array and a `log_lh: f64` field tracking the log-likelihood normalization constant. The constructor sets `log_lh: 0.0` with no parameter to override. Callers that need a non-zero initial `log_lh` must mutate the field after construction, violating the "initialize fully upfront" principle.

## Fix

Add `log_lh` as a constructor parameter, or provide a separate constructor that accepts it. Audit callers that mutate `log_lh` after construction and convert them to use the new parameter.

## Related issues

- Source: [N-dead-code-orphaned-apis.md](../issues/N-dead-code-orphaned-apis.md) -- delete after full resolution
