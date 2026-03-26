# Pseudo-count smoothing on initial mugration pi

v1 applies `apply_pseudo_counts(pi, pc)` to the equilibrium frequencies before constructing the initial mugration GTR model at [packages/treetime/src/commands/mugration/run.rs#L221](../../packages/treetime/src/commands/mugration/run.rs#L221). This adds `pc` (default 1.0) to each element and renormalizes, producing a smoother prior for the first reconstruction pass.

v0 passes raw weights (or uniform) directly to `GTR.custom(pi=weights, ...)` at [packages/legacy/treetime/treetime/wrappers.py#L766](../../packages/legacy/treetime/treetime/wrappers.py#L766) and reserves `pc` for `infer_gtr()` regularization only.

## Impact

For uniform pi (no weights), the pseudo-count addition has no effect (`[1/n + c, ...] / sum = [1/n, ...]`). For weight-based pi, the first reconstruction pass uses a flatter prior than v0, which shifts posterior probabilities at ambiguous internal nodes. The iterative refinement replaces the initial GTR with data-estimated parameters, limiting the impact to the first set of transition counts.

## Rationale

Pseudo-count smoothing prevents the initial model from having extreme equilibrium frequencies that could dominate the first posterior computation. This is standard Bayesian practice for categorical distributions (Dirichlet prior with concentration parameter `pc`).
