# Mixed-NaN distribution semantics are undefined

Sampled distributions can contain both finite values and `NaN`, but the distribution API does not define whether `NaN` means an invalid likelihood, a missing sample, or a value to ignore. Current reductions can fail inside `ndarray-stats`, and normalization then becomes empty.

## Evidence

- `fn neglog_function_normalize()` [packages/treetime-distribution/src/distribution_core/distribution.rs#L202-L207](../../packages/treetime-distribution/src/distribution_core/distribution.rs#L202-L207) converts a failed or non-finite minimum reduction into `Distribution::Empty`.

### NaN propagation through combined negative-log values

When combining negative-log amplitude arrays, `f64::min` ignores an isolated NaN operand (returns the non-NaN value). A NaN that survives reaches `neglog_function_normalize()`, where `ndarray-stats::QuantileExt::min()` returns `Err(MinMaxError::UndefinedOrder)` and the `.ok()` path normalizes the result to `Distribution::Empty`. Invalid numeric input is swallowed without any error or warning.

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
