# Timetree backward pass combines child time messages in plain probability space

## Summary

The timetree backward pass combines child node time-distribution messages with plain-probability multiplication and per-step peak normalization, not a log-scaled representation. Tail mass underflows to exactly zero before the coalescent factor is applied, so the coalescent contribution's log-space accumulation protects only its own dynamic range, not the product it multiplies into.

## Mechanism

In [`packages/treetime/src/timetree/inference/backward_pass.rs`](../../packages/treetime/src/timetree/inference/backward_pass.rs) child messages are folded with `distribution_multiplication(&current, parent_message)?.normalize()` under the `Plain` y-axis policy. Each step rescales the peak to 1 but lets the tails decay in plain `f64`; for deep trees or wide grids they reach exactly 0.

The coalescent factor is then applied by `distribution_apply_neg_log_weight` ([`packages/treetime-distribution/src/distribution_ops/log_cost.rs`](../../packages/treetime-distribution/src/distribution_ops/log_cost.rs)), which shifts by the combined minimum of `contribution(t) - ln(amplitude(t))`. Where the incoming amplitude has already underflowed to 0, the point carries no mass and is dropped via `Plain::is_defined`; the log-space shift cannot recover mass the plain-space product already lost.

This is a scientific-accuracy concern about lost tail mass, not a robustness crash. A separate crash on this path (a `~-1e-26` tail amplitude reaching `ln()` and producing `NaN`, collapsing the whole distribution to `Empty` under `--coalescent`) is fixed: `apply_neg_log_weight` now maps non-positive amplitudes to `+inf` cost instead of taking `ln`. That negative amplitude originates in the plain-space product's linear interpolation, not in this underflow ([N-distribution-function-product-negative-roundoff.md](N-distribution-function-product-negative-roundoff.md)). The bias described here remains open.

The crate provides `ScaledDistribution` (`log_scale * normalized inner`, [`packages/treetime-distribution/src/distribution_scaled/scaled.rs`](../../packages/treetime-distribution/src/distribution_scaled/scaled.rs)) for underflow-robust products, but the timetree backward pass does not use it.

## Impact

When the coalescent prior (or any factor) favors a node time in a region where the child-message product has underflowed, the posterior peak lands on residual near-zero mass, biasing or destroying the inferred node-time distribution.

This is distinct from the sequence-marginal asymmetry in [M-inference-forward-backward-asymmetry.md](M-inference-forward-backward-asymmetry.md): that entry concerns dense sequence profiles; this concerns node time distributions.

Occurrence is unverified. It requires wide $T_c$ or large lineage counts, where the coalescent contribution spans a large range over the message support. An empirical check -- driving such a case and dumping the pre- and post-coalescent distributions -- is needed to determine whether it triggers on real datasets.

## Related issues

- [M-distribution-product-grid-resolution-diverges-from-v0.md](M-distribution-product-grid-resolution-diverges-from-v0.md): the same product path resamples onto a uniform grid, a separate accuracy concern.
