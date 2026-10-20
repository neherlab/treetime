# GTR inference does not report non-convergence to its callers

`infer_gtr_from_counts()` returns `Ok(...)` with a warning log (`progress_warn!`) when the GTR inference does not converge [packages/treetime/src/gtr/infer_gtr/common.rs#L69-L79](../../packages/treetime/src/gtr/infer_gtr/common.rs#L69-L79). `struct InferGtrResult` [packages/treetime/src/gtr/infer_gtr/common.rs#L111](../../packages/treetime/src/gtr/infer_gtr/common.rs#L111) has no convergence flag, so callers cannot tell a converged model from a non-converged one without parsing the log.

## Fix

Add a convergence flag to `InferGtrResult` so that callers can detect and handle non-convergence, for example to report it in the run record.

## Validation

A test with counts that do not converge within the iteration limit reads the flag from the result.
