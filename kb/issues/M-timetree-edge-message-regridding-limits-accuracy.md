# Edge-message regridding limits the accuracy of the marginal likelihood

## Summary

The backward pass of time inference sends every edge message through `pub fn convolve_across_edge()`, which resamples the convolution onto a mass window: the window keeps the probability mass above `EPS = 5e-4` of each operand, and the result gets `GRID_POINTS = 300` points. Linear interpolation of the negative log onto these 300 points, and the soft tails that replace the cut mass, put an error of $10^{-6}$ to $10^{-5}$ into the log marginal likelihood of small test trees. The marginal clock rate step reads the rate from a parabola through three such likelihoods, so its target rate amplifies the error to about $10^{-3}$.

Locations:

- `pub(crate) const GRID_POINTS` and `pub(crate) const EPS` [packages/treetime/src/timetree/inference/runner.rs#L29-L36](../../packages/treetime/src/timetree/inference/runner.rs#L29-L36)
- `fn send_backward_message()` calls `convolve_across_edge` with them [packages/treetime/src/timetree/inference/backward_pass.rs#L160-L179](../../packages/treetime/src/timetree/inference/backward_pass.rs#L160-L179)
- `pub fn convolve_across_edge()` and its `resample_to_mass_window` call [packages/treetime-distribution/src/distribution_ops/edge_convolution.rs#L10-L43](../../packages/treetime-distribution/src/distribution_ops/edge_convolution.rs#L10-L43)

## Affected computation

The marginal clock rate probes: `pub(crate) fn step_marginal_clock_rate()` evaluates `fn marginal_log_lh_at_rate()` at the reference rate and at the rate times and divided by 1.25, and `fn log_rate_target()` places the target at the maximum of the parabola through the three values [packages/treetime/src/timetree/optimization/marginal_clock_rate.rs#L18-L115](../../packages/treetime/src/timetree/optimization/marginal_clock_rate.rs#L18-L115). The difference of nearby likelihoods divided by their curvature turns a likelihood error of $10^{-5}$ into a target rate error of $10^{-3}$.

## Measurements

All values are absolute errors against a closed-form or exactly integrated oracle, on fixtures fine enough that their own interpolation error is below the reported value.

| Test | Oracle | Error with `GRID_POINTS = 300`, `EPS = 5e-4` | Error with `GRID_POINTS = 900`, `EPS = 1e-9` |
|---|---|---|---|
| `test_marginal_clock_rate_matches_closed_form_curve`, rate 2, 2.5, 4 | $-\ln\mu - (\mu - 3)^2$ | 2.2e-5, 1.1e-6, 1.2e-5 | 1.1e-5 (largest) |
| `test_marginal_clock_rate_follow_slope_window_takes_a_maximum_past_the_probes` (target rate) | parabola through the closed form | 1.4e-3 | 3.1e-5 |
| `test_backward_pass_root_carries_the_marginal_likelihood` | analytic Gaussian integral | 8.3e-7 (fixture spacing 0.0005), 1.05e-6 (spacing 0.0001) | not measured |
| `test_marginal_likelihood_matches_piecewise_integration` | closed-form integration of each linear piece, converged by doubling to 3e-10 | 5.1e-7 and 8.7e-7 without a coalescent prior, 1.8e-6 to 3.3e-6 with one | not measured |

The two marginal clock rate tests carry the tolerances $10^{-4}$ and $10^{-2}$ with a message that names this issue.

> [!IMPORTANT]
> **Investigation required.** The gap of the marginal likelihood test with a coalescent prior (up to 3.3e-6) is larger than without one. Measure it with a finer edge grid to tell whether the regridding causes it, or the coalescent factor applied at the grid points of the subtree distribution, whose kinks at the prior node times fall between grid points.
>
> The tail laws of the coalescent factor are not the cause: `distribution_multiply_by_fn()` keeps the tail laws of its input ([M-distribution-multiply-by-fn-keeps-input-tail-laws.md](M-distribution-multiply-by-fn-keeps-input-tail-laws.md)), but a refit of the tails after the weight, and the slope of the weight added to the kept slope, both leave every gap of this test unchanged in all printed digits (largest 3.328e-6), because the grid edges of the multiplied distributions lie 36 to 69 negative-log units above their minimum.

## Fix direction

Decide the grid size and the mass window of the edge messages in the regrid round of the distribution design ([kb/proposals/distribution-log-space-and-hard-soft-boundaries-2.md](../proposals/distribution-log-space-and-hard-soft-boundaries-2.md), "Regridding"): a mass model whose cut mass enters the log marginal likelihood exactly, and a grid rule that sets the point count from the curvature of the negative log instead of a fixed count. Then tighten the tolerances above to at most 1e-6.
