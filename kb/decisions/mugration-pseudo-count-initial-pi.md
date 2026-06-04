# Pseudo-count smoothing on initial mugration pi

By default v1 matches v0: the initial mugration GTR model is built from the raw equilibrium frequencies, and the pseudo-count is reserved for `infer_gtr` regularization. v0 passes raw weights directly to `GTR.custom(pi=weights, ...)` at [packages/legacy/treetime/treetime/wrappers.py#L766](../../packages/legacy/treetime/treetime/wrappers.py#L766).

v1 offers an opt-in `--smooth-initial-pi` flag. When set, it applies `apply_pseudo_counts(pi, Some(pc.unwrap_or(1.0)))` to the equilibrium frequencies before constructing the initial GTR model at [packages/treetime/src/mugration/mugration.rs#L141-L143](../../packages/treetime/src/mugration/mugration.rs#L141-L143), producing a flatter prior for the first reconstruction pass. The flag uses the same effective pseudo-count as the refinement path, so the two pi-smoothing paths stay consistent.

`fixed_pi` (the weight-derived equilibrium pinned through refinement) is always captured from the raw frequencies at [packages/treetime/src/mugration/mugration.rs#L134](../../packages/treetime/src/mugration/mugration.rs#L134), before any smoothing. With smoothing off (the default) the initial pi equals `fixed_pi`, matching v0 exactly.

## Impact

For uniform pi (no weights) the flag is a no-op (`[1/n + c, ...] / sum = [1/n, ...]`). For weight-based pi, enabling the flag uses a flatter prior than v0 for the first reconstruction pass, shifting posterior probabilities at ambiguous internal nodes. The iterative refinement replaces the initial GTR with data-estimated parameters, limiting the effect to the first set of transition counts. The final equilibrium is always `fixed_pi`, independent of the flag.

## Rationale

Pseudo-count smoothing prevents the initial model from having extreme equilibrium frequencies that could dominate the first posterior computation. This is standard Bayesian practice for categorical distributions (Dirichlet prior with concentration parameter `pc`). It is opt-in because exact v0 parity is the default; the divergence occurs only when the flag is set.
