# Mixed-NaN distribution semantics are undefined

Sampled distributions can contain both finite values and `NaN`, but the distribution API does not define whether `NaN` means an invalid likelihood, a missing sample, or a value to ignore. Normalization and the most-likely-time search reject `NaN` with an error, while construction and the other reductions accept it, so the same distribution passes some operations and fails others.

## Evidence

- `fn neglog_function_normalize()` [packages/treetime-distribution/src/distribution_core/distribution.rs#L210](../../packages/treetime-distribution/src/distribution_core/distribution.rs#L210) reports an error when the minimum reduction meets `NaN`, and `fn neglog_peak()` [packages/treetime-distribution/src/distribution_core/distribution.rs#L228](../../packages/treetime-distribution/src/distribution_core/distribution.rs#L228) reports an error for a `NaN` or `-inf` peak. A peak of `+inf` (zero probability everywhere) normalizes to `Distribution::Empty`.
- `fn DistributionFunction.likely_time()` [packages/treetime-distribution/src/distribution_core/function.rs#L254](../../packages/treetime-distribution/src/distribution_core/function.rs#L254) and `fn DistributionFormula.likely_time()` [packages/treetime-distribution/src/distribution_core/formula.rs#L60](../../packages/treetime-distribution/src/distribution_core/formula.rs#L60) report an error when the values contain `NaN`.
- `DistributionFunction` construction accepts `NaN` ordinates, so a `NaN` sample is detected only when one of these operations runs, far from where it was produced.

### NaN propagation through combined negative-log values

When combining negative-log amplitude arrays, `f64::min` ignores an isolated NaN operand (returns the non-NaN value), so some `NaN` values disappear during arithmetic and others survive until normalization reports them.

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
