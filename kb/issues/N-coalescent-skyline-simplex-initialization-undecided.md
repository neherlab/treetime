# Skyline simplex initialization is scale-independent and one-sided

> **Resolved / obsolete.** `create_initial_simplex()` no longer exists. The skyline
> is now a convex solve in $z=\ln T_c$ by Newton's method, warm-started from the
> decoupled optimum $z_i=\ln(I_i/M_i)$ (no simplex). See
> [decisions/coalescent-skyline-convex-log-tc.md](../decisions/coalescent-skyline-convex-log-tc.md).

`fn create_initial_simplex()` adds $0.5$ to one log-$T_c$ coordinate per vertex regardless of grid spacing, parameter scale, or dimension, and never probes the negative direction [packages/treetime/src/coalescent/skyline.rs#L230-L248](../../packages/treetime/src/coalescent/skyline.rs#L230-L248). The value corresponds to multiplying one $T_c$ coordinate by $e^{0.5}$, but no convergence or conditioning contract justifies it.

## Design axes

### Scale

- O1. Derive each perturbation from the local parameter scale and optimizer tolerances.
- O2. Use a dimensionless fixed perturbation supported by a convergence study over the supported grid dimensions.

### Direction

- O1. Construct a balanced simplex around the initial estimate.
- O2. Keep a one-sided simplex only if boundary constraints require it and document the constraint.

## Recommendation

Use locally scaled, balanced vertices unless an independent convergence study demonstrates a better initialization. Create an implementation ticket after the scale and direction contracts are approved.

## Related issues

- [N-coalescent-skyline-quadrature-contract-undecided.md](N-coalescent-skyline-quadrature-contract-undecided.md)
