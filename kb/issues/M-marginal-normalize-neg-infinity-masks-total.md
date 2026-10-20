# `normalize_inplace` NEG_INFINITY masks all contributions

> [!IMPORTANT]
> **Decision required.** The finite-penalty, exclusion, and masked-accumulation options define different likelihoods. The intended probability model and convergence behavior must be derived before one is selected. See "Fix" and "Readiness".

## Summary

`normalize_inplace` returns `f64::NEG_INFINITY` for degenerate rows and substitutes uniform distributions. A single degenerate row makes the total log-likelihood negative infinity, masking contributions from all well-determined rows.

## Details

`fn normalize_inplace()` [packages/treetime/src/partition/marginal/shared/normalize.rs#L6-L21](../../packages/treetime/src/partition/marginal/shared/normalize.rs#L6-L21) iterates rows, normalizing each to sum to 1.0. For rows where the sum is zero or non-finite, it fills the row with a uniform distribution (1/n_cols) and adds `NEG_INFINITY` to the running `log_lh` total.

Since `log_lh` accumulates via addition, a single degenerate row produces:

```
log_lh = ... + finite_contributions + NEG_INFINITY = NEG_INFINITY
```

The callers in [packages/treetime/src/partition/marginal/shared/pass.rs](../../packages/treetime/src/partition/marginal/shared/pass.rs) treat the result differently:

- Forward message to a child ([pass.rs#L298-L300](../../packages/treetime/src/partition/marginal/shared/pass.rs#L298-L300)): `fn forward_log_lh_add_normalization()` [packages/treetime/src/partition/marginal/shared/normalize.rs#L37-L43](../../packages/treetime/src/partition/marginal/shared/normalize.rs#L37-L43) drops a `NEG_INFINITY` normalization, so the message keeps a finite log-likelihood
- Root profile ([pass.rs#L189-L192](../../packages/treetime/src/partition/marginal/shared/pass.rs#L189-L192)) and dense node profile ([pass.rs#L305-L309](../../packages/treetime/src/partition/marginal/shared/pass.rs#L305-L309)): the result is added directly, so one degenerate row makes the node log-likelihood `NEG_INFINITY`

Convergence checks downstream see `NEG_INFINITY` even when the rest of the tree is well-determined.

A row is degenerate only when its unnormalized profile is zero or non-finite. A gap or unknown character does not cause this: the alphabet maps both to an all-ones profile ([packages/treetime/src/alphabet/alphabet_config.rs#L45](../../packages/treetime/src/alphabet/alphabet_config.rs#L45), [#L56](../../packages/treetime/src/alphabet/alphabet_config.rs#L56)), so an all-gap column normalizes to a finite value. The passes clamp branch lengths to a positive minimum, so propagated messages are strictly positive. A degenerate row then needs an underflow or a non-finite input value. No reproduction on real input is known.

## Impact

- Convergence criteria based on likelihood changes become unreliable when a degenerate row occurs
- A single degenerate position masks the likelihood signal of the whole tree
- Optimizer cannot distinguish "improving" from "stuck" when total is negative infinity

## Fix

Options:

1. Return a finite penalty for degenerate rows (e.g., `ln(1/n_cols) * n_positions`) instead of `NEG_INFINITY`
2. Track degenerate row count separately and exclude from convergence arithmetic
3. Use a masked accumulator that skips degenerate contributions

## Readiness

The listed finite-penalty, exclusion, and masked-accumulation choices define different likelihoods. Derive the intended probability model and convergence behavior before selecting a correction.
