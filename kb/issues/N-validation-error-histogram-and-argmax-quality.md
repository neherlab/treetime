# Validation crate histogram and argmax quality issues

Two quality issues in the validation crate's metrics code.

### Error histogram functions share 80% logic

`compute_error_histogram` and `compute_error_histogram_signed` differ only in filter predicate (`x > 0.0` vs `x.is_finite()`) and the `use_log_scale` parameter (unsigned applies `log10()` transform; signed uses raw values). The remaining 25 lines of binning, edge construction, center computation, and return structure are identical.

### argmax uses unspecified tie behavior

Two call sites use `ndarray_stats::argmax()` which has unspecified tie-breaking behavior, instead of the project's `argmax_first()` which guarantees deterministic first-index tie-breaking for Python parity.

## Locations

- `packages/treetime-validation/src/testing/metrics/distribution/histograms.rs:75-133:` `compute_error_histogram()`
- `packages/treetime-validation/src/testing/metrics/distribution/histograms.rs:135-177:` `compute_error_histogram_signed()`
- `packages/treetime-validation/src/testing/metrics/aggregate/domain_agreement/peak_metrics.rs:27-28:` two `argmax()` calls

## Impact

Histogram duplication: maintainability. Argmax: on continuous analytical distributions ties are vanishingly improbable, so this is not a correctness risk in practice.

## Action

Unify histogram functions with a config parameter (signed/unsigned mode). Replace `argmax()` with `argmax_first()`.
