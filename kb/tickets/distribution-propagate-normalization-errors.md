# Propagate distribution normalization errors

Separate explicit domain emptiness from formula-evaluation and grid-construction failures during negative-log normalization.

## Required changes

- Return `Result<Distribution<Plain>, Report>` from `Distribution<NegLog>::to_plain_normalized()` [packages/treetime-distribution/src/distribution_core/distribution.rs#L378](../../packages/treetime-distribution/src/distribution_core/distribution.rs#L378).
- Return `Ok(Distribution::Empty)` only for the explicit `Distribution::Empty` variant.
- Make `neglog_function_to_plain_normalized()` [packages/treetime-distribution/src/distribution_core/distribution.rs#L437](../../packages/treetime-distribution/src/distribution_core/distribution.rs#L437) fallible and propagate `DistributionFunction::from_start_dx_values()` errors with grid context.
- Propagate formula discretization errors with formula bounds and evaluation context.
- Update all callers, including [packages/treetime/src/timetree/inference/backward_pass.rs#L58](../../packages/treetime/src/timetree/inference/backward_pass.rs#L58), to propagate normalization failure instead of storing an empty distribution.

## Validation

- Formula-evaluation error, grid-construction error, and explicit-empty cases.
- Property tests for maximum one, non-negativity, finiteness, representable likelihood-ratio preservation, and common-offset invariance.
- A timetree backward-pass test proving a normalization error reaches the caller and no node distribution is committed.
- Existing underflow golden masters and full lint/test suite.

## Related issues

- Source: [kb/issues/M-distribution-normalization-erases-errors.md](../issues/M-distribution-normalization-erases-errors.md)
