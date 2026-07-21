# Mixed-NaN distribution semantics are undefined

Sampled distributions can contain both finite values and `NaN`, but the distribution API does not define whether `NaN` means an invalid likelihood, a missing sample, or a value to ignore. Current reductions produce operation-dependent results: `ndarray-stats` reductions can fail, normalization can become empty, and maximum queries can fall back to zero.

## Evidence

- `Distribution<Plain>::max_value()` [packages/treetime-distribution/src/distribution_core/distribution.rs#L193](../../packages/treetime-distribution/src/distribution_core/distribution.rs#L193) converts reduction and formula-discretization errors to `0.0`, conflating invalid input with a legitimate zero maximum.
- `fn neglog_function_to_plain_normalized()` [packages/treetime-distribution/src/distribution_core/distribution.rs#L437](../../packages/treetime-distribution/src/distribution_core/distribution.rs#L437) converts a failed minimum reduction into `Distribution::Empty`.
- `Distribution<Plain>::normalize()` [packages/treetime-distribution/src/distribution_core/distribution.rs#L213](../../packages/treetime-distribution/src/distribution_core/distribution.rs#L213) treats a non-finite maximum as empty.

### NaN propagation through combined negative-log values

When combining negative-log amplitude arrays, `f64::min` ignores an isolated NaN operand (returns the non-NaN value). The surviving NaN propagates through exponentiation into the plain-space array. `ndarray-stats::QuantileExt::max()` then returns `Err(MinMaxError::UndefinedOrder)`, and the existing `.ok().unwrap_or(0.0)` path normalizes the result to `Distribution::Empty`. Invalid numeric input is swallowed without any error or warning.

## Decision axes

### Meaning of `NaN` in sampled values

- **Reject any `NaN`:** preserves strict likelihood semantics and prevents partial evidence from being silently discarded.
- **Ignore `NaN` when finite samples remain:** supports missing samples, but an ordinary `Array1<f64>` cannot distinguish intentional missingness from numerical failure.
- **Represent missing samples explicitly:** use a mask or typed sample state and reject unmarked `NaN`. This supports incomplete grids without overloading floating-point error values.

Recommendation: reject every `NaN` in the current untyped arrays. Use an explicit representation if missing samples are a required domain feature.

### Enforcement boundary

- **Validate during `DistributionFunction` construction:** invalid states cannot enter the distribution API, and every downstream reduction shares one policy.
- **Validate in each reduction and conversion:** different operations can adopt different semantics, but the same distribution may be accepted by one operation and rejected by another.

Recommendation: validate during construction. Constructors are the boundary where unstructured samples become a valid distribution.

## Recommendation

Reject every `NaN` during construction unless the application requires missing samples; in that case, represent missingness explicitly and still enforce the contract during construction. This issue remains ticketless until the missing-sample requirement and representation are approved.

## Related issues

- [M-distribution-normalization-erases-errors.md](M-distribution-normalization-erases-errors.md)
