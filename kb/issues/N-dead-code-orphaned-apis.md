# Dead code and orphaned APIs

## Summary

Constructor hard-codes a field that callers must mutate post-construction.

## Instances

### DenseSeqDis::new() hard-codes log_lh: 0.0

`packages/treetime/src/representation/payload/dense.rs:53-55:`

`struct DenseSeqDis` holds a 2D probability array and a `log_lh: f64` field tracking the log-likelihood normalization constant. The constructor sets `log_lh: 0.0` with no parameter to override. Callers that need a non-zero initial `log_lh` must mutate the field after construction, violating the "initialize fully upfront" principle.

## Impact

- Constructor requiring post-construction mutation is error-prone
