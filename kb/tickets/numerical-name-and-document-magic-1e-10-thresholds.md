# Name and document magic 1e-10 thresholds in timetree and relaxed_clock

Two locations use hardcoded `1e-10` thresholds without named constants, documentation, or validation of degenerate inputs.

## Instances

### create_poisson_branch_distributions hardcoded threshold

v1: [`packages/treetime/src/commands/timetree/utils.rs#L100`](../../packages/treetime/src/commands/timetree/utils.rs#L100)

Hardcoded `1e-10` threshold with no `n_points==1` or `mu<=0` validation. Mishandles `branch_length==0` (produces degenerate grid).

### 1e-10 magic denominator in relaxed_clock.rs

v1: [`packages/treetime/src/commands/timetree/optimization/relaxed_clock.rs#L70`](../../packages/treetime/src/commands/timetree/optimization/relaxed_clock.rs#L70), [`relaxed_clock.rs#L91`](../../packages/treetime/src/commands/timetree/optimization/relaxed_clock.rs#L91), [`relaxed_clock.rs#L106`](../../packages/treetime/src/commands/timetree/optimization/relaxed_clock.rs#L106)

Three locations use `1e-10` as a denominator floor. No named constant. Related to the unit mismatch in `H-timetree-relaxed-clock-unit-mismatch.md`.

## Impact

- Degenerate inputs (zero branch length, mu<=0) not rejected at boundaries
- Magic constants without named constants reduce maintainability

## Related issues

- Source: [N-numerical-stability-magic-constants.md](../issues/N-numerical-stability-magic-constants.md) -- delete after full resolution
