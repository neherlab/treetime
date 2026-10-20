# Multiplying a distribution by a function keeps the tail laws of the input

## Summary

`pub fn distribution_multiply_by_fn()` adds the negative log of a weight $w(t)$ to the grid ordinates of a negative-log distribution $f(t)$, but copies the left and right tail laws of $f$ unchanged. Off the grid, the product then follows the tail of $f$ instead of the tail of $f + w$. The tail of a soft side is the linear law $y(t) = y_\text{edge} + s\,(t - t_\text{edge})$, whose mass is $e^{-y_\text{edge}} / |s|$, so a kept slope $s_f$ in place of $s_f + s_w$ misstates the tail mass by the factor $|s_f + s_w| / |s_f|$.

Example: the input $f(t) = -2t$ with a soft left tail and the weight $w(t) = 0.5t$ give the product $-1.5t$. One unit past the left edge, at $t = -1$, the kept tail gives $2$ instead of $1.5$.

Locations:

- `pub fn distribution_multiply_by_fn()` [packages/treetime-distribution/src/distribution_ops/multiply_by_fn.rs#L8-L48](../../packages/treetime-distribution/src/distribution_ops/multiply_by_fn.rs#L8-L48); the tail laws of the input are copied at [#L39-L41](../../packages/treetime-distribution/src/distribution_ops/multiply_by_fn.rs#L39-L41)
- The test `test_multiply_by_fn_keeps_the_incoming_tail_laws` in `packages/treetime-distribution/src/distribution_ops/__tests__/test_multiply_by_fn.rs` asserts the kept slope $-2$ for this example

The function refitted both soft tails from the outermost grid points after it added the weight before commit `0953dfadd` ("keep the incoming tails when multiplying by a function"). The coalescent cost has a step where the number of lineages changes. A step inside the fit points tilted the fitted left tail outward, the clamp of `SoftTailLaw::fit()` made it flat, and the node distribution got infinite mass (3 to 7 recent nodes per simulated run with a coalescent prior). Keeping the incoming tail removed these failures. A plain refit is therefore not the fix.

## Consumers

- The coalescent factor of internal nodes, `fn apply_coalescent_prior()` [packages/treetime/src/timetree/inference/backward_pass.rs#L121-L134](../../packages/treetime/src/timetree/inference/backward_pass.rs#L121-L134)
- The coalescent factor of tips, `fn apply_leaf_coalescent_factor()` [packages/treetime/src/timetree/inference/backward_pass.rs#L150-L173](../../packages/treetime/src/timetree/inference/backward_pass.rs#L150-L173)
- Through them: the tail mass in `pub fn Distribution::log_mass()` [packages/treetime-distribution/src/distribution_core/distribution.rs#L216](../../packages/treetime-distribution/src/distribution_core/distribution.rs#L216), the log marginal likelihood of node times, and the marginal clock rate step `pub(crate) fn step_marginal_clock_rate()` [packages/treetime/src/timetree/optimization/marginal_clock_rate.rs#L18](../../packages/treetime/src/timetree/optimization/marginal_clock_rate.rs#L18), which reads the rate from a parabola through three such likelihoods

## Measurement

A temporary change of `distribution_multiply_by_fn()` was measured against `test_marginal_likelihood_matches_piecewise_integration` (`packages/treetime/src/timetree/inference/__tests__/test_marginal_likelihood.rs`), whose oracle integrates each linear piece in closed form. The gap is the log marginal likelihood of the backward pass minus the oracle value:

| Case | Kept tails (current) | Refit after the weight | Slope of the weight added to the kept slope |
|---|---|---|---|
| rate 1e-3, no prior | 5.059e-7 | 5.059e-7 | 5.059e-7 |
| rate 5e-4, no prior | -8.716e-7 | -8.716e-7 | -8.716e-7 |
| rate 1e-3, $T_c = 5$ | 2.222e-6 | 2.222e-6 | 2.222e-6 |
| rate 5e-4, $T_c = 5$ | 3.328e-6 | 3.328e-6 | 3.328e-6 |
| rate 2e-3, $T_c = 20$ | 1.774e-6 | 1.774e-6 | 1.774e-6 |

The composed variant changed the left slope as expected: in the case rate 5e-4 with $T_c = 5$, from $-0.869$ to $-1.069$ and from $-1.299$ to $-1.499$, the slope $-1/T_c$ of the weight. The gaps did not change in any digit, because the grid edges of the multiplied distributions lie 36 to 69 negative-log units above their minimum, so the tail mass is below $e^{-36} \approx 2 \times 10^{-16}$. The tail law is therefore not the cause of the larger gap with a coalescent prior in this fixture (see [M-timetree-edge-message-regridding-limits-accuracy.md](M-timetree-edge-message-regridding-limits-accuracy.md)). The defect matters where a grid edge lies close to the peak, for example after the mass window of an edge message cuts at `EPS = 5e-4` of the mass, where $y_\text{edge}$ is about $\ln(1/\text{EPS}) \approx 7.6$ above the minimum.

> [!IMPORTANT]
> **Investigation required.** No test measures the effect of the kept tail law on a distribution whose grid edge is near its peak. Build a fixture where the multiplied distribution has a soft tail with mass above 1e-6, and compare `log_mass()` with the closed-form mass of $f + w$.

## Fix direction

Compose the weight with the tail laws instead of refitting from the grid: on each soft side, add the slope of the weight beyond the grid edge to the slope of the incoming law (`SoftTailLaw::compose_multiply()`), and keep the law decaying. The slope of the weight must come from the weight past the edge, not from the grid points, so a step inside the fit points does not tilt the tail. Do this inside the mass model of the next regrid round of the distribution design ([kb/proposals/distribution-log-space-and-hard-soft-boundaries-2.md](../proposals/distribution-log-space-and-hard-soft-boundaries-2.md), "Regridding"), which decides how cut mass enters the log marginal likelihood. Then replace the assertion of the kept slope in `test_multiply_by_fn_keeps_the_incoming_tail_laws` with the slope of the product, $-1.5$ in the example.
