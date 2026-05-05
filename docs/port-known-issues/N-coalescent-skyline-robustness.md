# Coalescent skyline robustness deficiencies

## Summary

Seven defects in the coalescent skyline subsystem affecting numerical robustness: silent skipping of timeless nodes, unclamped optimizer output, ignored convergence flags, non-scaling simplex perturbation, flat extrapolation outside grid, first-order integration approximation, and sign-convention inconsistency.

## Instances

### collect_tree_events silently skips nodes without likely_time

`packages/treetime/src/commands/timetree/coalescent/events.rs:28-41:`

Lineage count k(t) is wrong for nodes with `time_distribution` but `None` for `likely_time`. These nodes are silently excluded from the event timeline, under-counting active lineages.

### Skyline unclamped return from optimizer

`packages/treetime/src/commands/timetree/coalescent/skyline.rs:154:`

`log_tc_values` returns raw optimizer output while the cost function internally clamps to [-200, 100]. The optimizer can return values outside the clamp range that were never evaluated by the actual cost function.

### optimize_tc ignores argmin convergence flag

`packages/treetime/src/commands/timetree/coalescent/optimize_tc.rs:75:`

Reports `success=true` when Brent's method returns `None` (no convergence). Callers cannot distinguish converged from unconverged results.

### Simplex perturbation fixed at +0.5

`packages/treetime/src/commands/timetree/coalescent/skyline.rs:310:`

`fn create_initial_simplex` perturbs each dimension of the initial Nelder-Mead vertex by a fixed +0.5 in log-Tc space (~1.65x multiplicative factor). The perturbation does not scale with the number of grid points or their spacing, and is asymmetric (positive direction only).

### Tc distribution flat extrapolation outside grid

`packages/treetime/src/commands/timetree/coalescent/skyline.rs:233-238:`

Calls to `tc_dist.eval(t)` at times outside the lineage-count event range return the boundary value silently, implying constant Tc beyond observed data without any indication to callers.

### Midpoint integration for non-constant Tc

`packages/treetime/src/commands/timetree/coalescent/integration.rs:63-68:`

First-order midpoint approximation for the coalescent integral. Error up to 2.5% for rapidly varying Tc between grid points. Higher-order quadrature (Simpson's or Gauss-Legendre) would reduce this.

### Leaf contribution sign-convention inconsistency

`packages/treetime/src/commands/timetree/coalescent/contributions.rs:89-105:`

Return type is `DistributionNegLog` but values follow log-probability convention (positive = more likely). The type name implies negated log-probability (positive = less likely). Related to the time notation conflict documented in `M-coalescent-time-notation-conflict.md`.

## Impact

- Lineage counts underestimated when nodes lack `likely_time`
- Optimizer results unreliable when convergence not checked
- Skyline Tc estimates biased by non-scaling simplex and flat extrapolation
- Integration error accumulates for rapidly changing population sizes
