# Numerical stability magic constants

## Summary

Six locations use hardcoded numeric thresholds without named constants, documentation, or validation of degenerate inputs.

## Instances

### create_poisson_branch_distributions hardcoded threshold

`packages/treetime/src/commands/timetree/utils.rs:100:`

Hardcoded `1e-10` threshold with no `n_points==1` or `mu<=0` validation. Mishandles `branch_length==0` (produces degenerate grid).

### 1e-10 magic denominator in relaxed_clock.rs

`packages/treetime/src/commands/timetree/optimization/relaxed_clock.rs:70,91,106:`

Three locations use `1e-10` as a denominator floor. No named constant. Related to the unit mismatch in `H-timetree-relaxed-clock-unit-mismatch.md`.

### Discrete MIN_POSITIVE clamp biases near-zero divisor

`packages/treetime/src/representation/partition/discrete.rs:207:`

Silently biases near-zero values to `f64::MIN_POSITIVE` without propagating information about the degenerate log-likelihood to callers.

### SUPERTINY_NUMBER distorts column-stochastic normalization

`packages/treetime/src/gtr/infer_gtr/dense.rs:199:`

`SUPERTINY_NUMBER` (1e-24) added to `expQt` result distorts column-stochastic normalization. The additive constant prevents exact zeros but shifts column sums away from 1.0.

### Non-positive Tc yields NaN in contributions

`packages/treetime/src/commands/timetree/coalescent/contributions.rs:117-178:`

`fn compute_internal_contribution_single()` has no guard on non-positive Tc. When Tc(t) <= 0, lambda(t) goes negative or infinite, and ln() yields NaN or negative infinity. Callers do not validate Tc positivity before invoking.

### debug_assert!(gamma > 0.0) not a production guard

`packages/treetime/src/commands/timetree/inference/branch_length_likelihood.rs:61:`

The assertion is stripped in release builds. When gamma <= 0 reaches this code path in production, the computation produces inf/NaN silently.

## Impact

- Silent NaN/inf propagation in production builds for edge-case inputs
- Column-stochastic property violated by additive perturbation
- Degenerate inputs (zero branch length, non-positive Tc, zero gamma) not rejected at boundaries
