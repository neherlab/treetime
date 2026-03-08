# Coalescent likelihood always returns None

`compute_coalescent_likelihood()` is a stub that always returns `None`.
Coalescent likelihood is never computed, making convergence metrics incomplete
when using coalescent models.

## Location

[`likelihood.rs#L77`](../../packages/treetime/src/commands/timetree/convergence/likelihood.rs#L77)

## Impact

- `lh_coal` is always `None` in convergence metrics
- `lh_total` excludes coalescent contributions
- Convergence detection cannot account for coalescent model fit

## Related issues

- [Positional likelihood metric differs from v0](M-timetree-positional-likelihood-metric.md)
  another convergence metric gap
- [Internal node dates missing for coalescent modes](M-timetree-internal-dates-missing-coalescent.md)
  tracelog shows `lh_pos=NaN` for coalescent runs
