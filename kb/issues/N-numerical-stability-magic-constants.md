# Numerical stability magic constants

## Summary

Six locations use hardcoded numeric thresholds without named constants, documentation, or validation of degenerate inputs.

## Instances

### create_poisson_branch_distributions hardcoded threshold

`packages/treetime/src/timetree/utils.rs:100:`

Hardcoded `1e-10` threshold with no `n_points==1` or `mu<=0` validation. Mishandles `branch_length==0` (produces degenerate grid).

### 1e-10 magic denominator in relaxed_clock.rs

`packages/treetime/src/timetree/optimization/relaxed_clock.rs:70,91,106:`

Three locations use `1e-10` as a denominator floor. No named constant.

### Shared marginal forward pass MIN_POSITIVE clamp biases near-zero divisor

`packages/treetime/src/partition/marginal_core.rs:163:`

Silently biases near-zero values to `f64::MIN_POSITIVE` without propagating information about the degenerate log-likelihood to callers. Used by both dense and discrete partitions.

### SUPERTINY_NUMBER distorts column-stochastic normalization

`packages/treetime/src/gtr/infer_gtr/common.rs:199:`

`SUPERTINY_NUMBER` (1e-24) added to `expQt` result distorts column-stochastic normalization. The additive constant prevents exact zeros but shifts column sums away from 1.0.

### Non-positive Tc yields NaN in contributions

`packages/treetime/src/coalescent/contributions.rs:117-178:`

`fn compute_internal_contribution_single()` has no guard on non-positive Tc. When Tc(t) <= 0, lambda(t) goes negative or infinite, and ln() yields NaN or negative infinity. Callers do not validate Tc positivity before invoking.

### debug_assert!(gamma > 0.0) not a production guard

`packages/treetime/src/timetree/inference/branch_length_likelihood.rs:61:`

The assertion is stripped in release builds. When gamma <= 0 reaches this code path in production, the computation produces inf/NaN silently.

## Impact

- Silent NaN/inf propagation in production builds for edge-case inputs
- Column-stochastic property violated by additive perturbation
- Degenerate inputs (zero branch length, non-positive Tc, zero gamma) not rejected at boundaries

## Related tickets

- [kb/tickets/numerical-guard-non-positive-tc-in-coalescent-contributions.md](../tickets/numerical-guard-non-positive-tc-in-coalescent-contributions.md)
- [kb/tickets/numerical-min-positive-clamp-biases-near-zero-divisor.md](../tickets/numerical-min-positive-clamp-biases-near-zero-divisor.md)
- [kb/tickets/numerical-name-and-document-magic-1e-10-thresholds.md](../tickets/numerical-name-and-document-magic-1e-10-thresholds.md)
- [kb/tickets/numerical-promote-debug-assert-to-production-guards.md](../tickets/numerical-promote-debug-assert-to-production-guards.md)
- [kb/tickets/numerical-supertiny-number-distorts-column-stochastic-normalization.md](../tickets/numerical-supertiny-number-distorts-column-stochastic-normalization.md)
