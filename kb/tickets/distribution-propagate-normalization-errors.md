# Propagate distribution normalization errors

Separate explicit domain emptiness from formula-evaluation and grid-construction failures during negative-log normalization.

## Required changes

- Return `Result<Distribution<NegLog>, Report>` from `Distribution<NegLog>::normalize()` [packages/treetime-distribution/src/distribution_core/distribution.rs#L177-L199](../../packages/treetime-distribution/src/distribution_core/distribution.rs#L177-L199).
- Return `Ok(Distribution::Empty)` only for the explicit `Distribution::Empty` variant.
- Make `neglog_function_normalize()` [packages/treetime-distribution/src/distribution_core/distribution.rs#L202-L207](../../packages/treetime-distribution/src/distribution_core/distribution.rs#L202-L207) fallible and report a failed or non-finite minimum with grid context.
- Propagate formula discretization errors with formula bounds and evaluation context.
- Update all callers, including [packages/treetime/src/timetree/inference/backward_pass.rs#L41-L44](../../packages/treetime/src/timetree/inference/backward_pass.rs#L41-L44), to propagate normalization failure instead of storing an empty distribution.

## Validation

- Formula-evaluation error, grid-construction error, and explicit-empty cases.
- Property tests for minimum zero, non-negativity, finiteness, difference preservation, and common-offset invariance.
- A timetree backward-pass test proving a normalization error reaches the caller and no node distribution is committed.
- Existing underflow golden masters and full lint/test suite.

## Related issues

- Source: [kb/issues/M-distribution-normalization-erases-errors.md](../issues/M-distribution-normalization-erases-errors.md)
