# GTR inference does not report non-convergence to its callers

`fn infer_gtr_group()`, which `fn infer_gtr_from_counts()` wraps, returns `Ok(...)` with a warning log (`progress_warn!`) when the GTR inference does not converge or does not reach normalization [packages/treetime/src/gtr/infer_gtr/common.rs#L133-L143](../../packages/treetime/src/gtr/infer_gtr/common.rs#L133-L143). `struct InferGtrResult` [packages/treetime/src/gtr/infer_gtr/common.rs#L231](../../packages/treetime/src/gtr/infer_gtr/common.rs#L231) has no convergence flag, so callers cannot tell a converged model from a non-converged one without parsing the log.

## Fix

Add a convergence flag to `InferGtrResult` so that callers can detect and handle non-convergence, for example to report it in the run record.

## Validation

A test with counts that do not converge within the iteration limit reads the flag from the result.

## Proposed implementation plan

The warnings are in `fn infer_gtr_group()` at the cited lines; `struct InferGtrResult` holds only `W`, `pi`, `mu`. Users: `packages/treetime/src/partition/fitch/gtr_inference.rs`, `packages/treetime/src/gtr/refinement.rs`, and tests in `packages/treetime/src/gtr/infer_gtr/__tests__/`.

1. Add an enum `GtrConvergence { Converged, NotConverged, NotNormalized }` (or a struct with the final `pi` distance and iteration count) to `InferGtrResult` in `packages/treetime/src/gtr/infer_gtr/common.rs`, set from the same conditions as the two warnings. One status applies to all models of a group, because the check uses the maximum distance over the group.
2. Keep the warning log, and pass the status up: `fn gtr_from_inference()` and `refinement.rs` carry it to the run record or the model output where the command reports model fit details (check sibling commands for where fit diagnostics go).
3. Update every construction site of `InferGtrResult` (tests included).

Tests: in `test_group.rs`, one rstest case with counts that converge (status `Converged`) and one with `max_iter = 1` in `InferGtrOptions` on counts far from equilibrium (status `NotConverged`). Oracle: the convergence condition `max_distance(pi_old, pi) > dp` evaluated by hand on the chosen counts.

Validation: `just lint-rs`, `just test-rs`; `just generated-check` if the status enters a schema.
