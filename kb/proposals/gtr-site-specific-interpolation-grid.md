# Proposal: max-rate interpolation grid for site-specific GTR

## Summary

Evaluate an error-controlled alternative to the mean-rate interpolation grid used by site-specific GTR models. Replacing the mean by the maximum expected rate is a candidate heuristic, not an error bound.

## Motivation

The interpolation grid spacing controls the accuracy of the linear interpolation approximation. Grid spacing is set from the mean substitution rate across all sites. For heterogeneous alignments where a few sites evolve much faster than the mean, the grid is too coarse for those fast sites.

For the effective site generator $G_a=\mu_aQ_a$, the linear-interpolation remainder satisfies a bound of the form

$$\left\|P_a(t)-I_hP_a(t)\right\|\leq\frac{h^2}{8}\sup_\xi\left\|G_a^2e^{G_a\xi}\right\|.$$

The expected substitution rate can correlate with this error, but neither its mean nor its maximum generally bounds $\lVert G_a\rVert$, $\lVert G_a^2\rVert$, or the largest eigenvalue magnitude. A spectral scale such as $\max_{a,j}|\mu_a\lambda_{ja}|$ is more directly connected to curvature, although even that does not establish an entrywise probability-error target.

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

One candidate is to use `max(average_rate)` instead of `mean(average_rate)`. This contracts the covered branch-length interval and narrows its grid cells, but does not by itself ensure an accuracy target. A production change requires comparison with direct matrix exponentials over the supported generator and branch-length domain.

## Tradeoffs

- The retained 61-point array has unchanged memory use
- A larger scale contracts the physical interpolation interval and triggers exact fallback at shorter branch lengths
- Accuracy and runtime effects require separate measurement
- Diverges from v0 behavior
- For uniform-rate models (all sites equal), no change

## Decision needed

The mean-rate grid defines v0 parity. Current tests compare interpolated and direct transition matrices at `1e-2`; approximate column sums use `1e-8`, non-negativity uses `1e-10`, and equilibrium preservation uses `1e-2`. The golden master validates parity with v0's interpolation construction, not absolute interpolation accuracy. Changing the grid requires an approved numerical error contract and evidence across the supported domain.
