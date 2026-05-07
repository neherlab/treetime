# Proposal: max-rate interpolation grid for site-specific GTR

## Summary

Use the maximum per-site rate instead of the mean rate to calibrate the interpolation grid for site-specific GTR models.

## Motivation

The interpolation grid spacing controls the accuracy of the linear interpolation approximation. Grid spacing is set from the mean substitution rate across all sites. For heterogeneous alignments where a few sites evolve much faster than the mean, the grid is too coarse for those fast sites.

Linear interpolation error scales with `h^2 * |d^2/dt^2 exp(Q_a * t)|`, where `h` is the grid spacing. The second derivative depends on `Q_a^2`, which is controlled by the per-site rate and eigenvalue magnitude. Fast sites have larger second derivatives and need denser grids.

## Current behavior (v0 parity)

Both v0 and v1 use `mean(average_rate)` to set `rate_scale`:

```python
# v0: gtr_site_specific.py:335
self.rate_scale = self.average_rate().mean()
```

```rust
// v1: gtr_site_specific.rs
let rate_scale = (avg_rates.sum() / self.seq_len as f64).max(1e-10);
```

The grid extends to `10 / rate_scale` with 61 points. The fallback to exact computation activates when `t * rate_scale >= 10`.

## Proposed change

Use `max(average_rate)` instead of `mean(average_rate)` for grid calibration. This ensures the fastest site is within the accuracy target. Slower sites would have a denser-than-necessary grid (wastes some memory but not accuracy).

## Tradeoffs

- Denser grid for most sites (slight memory increase, no accuracy cost)
- Faster sites get accurate interpolation instead of relying on the fallback
- Diverges from v0 behavior
- For uniform-rate models (all sites equal), no change

## Decision needed

v0 has shipped with mean-rate grid for years. No evidence of practical accuracy problems. The proposal improves worst-case accuracy at the cost of v0 divergence. The property test validates interpolation invariants (stochastic, non-negative, equilibrium) at 1e-8 regardless of grid choice.
