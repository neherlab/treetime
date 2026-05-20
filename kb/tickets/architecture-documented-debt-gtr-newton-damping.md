# Architectural debt (documented, low priority)

## Summary

Three documented architectural choices that are suboptimal but stable: mixed stochastic conventions in GTR, looser-than-typical Newton tolerance, and a hardcoded damping floor.

## Instances

### Mixed row/column stochastic convention in GTR

`packages/treetime/src/gtr/gtr.rs:47-54:`

The `Q` matrix and `expQt` use different stochastic conventions in different contexts. Documented in code comments. The convention difference is intentional for performance (avoiding transpositions) but increases cognitive load for contributors working on GTR-related code.

### NEWTON_REL_TOL=0.001 looser than typical

`packages/treetime/src/optimize/method_newton.rs:10:`

Relative tolerance of 0.001 for Newton's method convergence. Documented and tested. Looser than the typical 1e-6 to 1e-8 for Newton's method. The choice is deliberate: branch length optimization does not require machine-precision convergence, and the looser tolerance reduces iteration count.

### apply_damping hardcodes 1% floor

`packages/treetime/src/commands/optimize/run.rs:501:`

`DAMPING_FLOOR=0.01`. Prevents damping from suppressing updates entirely. Documented as a named constant. The value matches v0 behavior and is not configurable.

## Impact

- Low: all three are documented, tested, and stable
- Mixed stochastic convention increases onboarding cost for GTR contributors
- Newton tolerance and damping floor values are not configurable but satisfy current use cases

## Related issues

- Source: [N-architectural-debt-documented.md](../issues/N-architectural-debt-documented.md) -- delete after full resolution
