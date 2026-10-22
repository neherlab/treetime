# Low-signal clock rate bias from first-order convolution error at grid edges

Date: 2026-10-21

> [!WARNING]
> **Needs review.** The convergence claim is reproduced: the rectangle sum converges at first order at hard grid ends and the trapezoid rule at second order ([kb/decisions/distribution-convolution-trapezoid-rule.md](../decisions/distribution-convolution-trapezoid-rule.md)). The numbers on 400 tips (the root log mass errors, the maximum of $\log Z$ at a quarter of the true rate, the MCMC maximum, and the final rates of the rate step) are not reproduced, and the MCMC used for them had a biased root move: [kb/issues/M-timetree-low-signal-rate-bias-claims-not-reproduced.md](../issues/M-timetree-low-signal-rate-bias-claims-not-reproduced.md).

## Problem

At low temporal signal, the log marginal likelihood of the node times (`log_lh_time`) under a coalescent prior rose toward slow clock rates. On validation `const` sets with $L=10^3$, $\mu=3\cdot10^{-4}$ and $T_c=10$, its maximum lay at about a quarter of the true rate in 10 of 10 replicates. A rate that followed the slope every round ran down to about $8\cdot10^{-5}$. The rate window of [timetree-clock-rate-by-marginal-likelihood.md](../decisions/timetree-clock-rate-by-marginal-likelihood.md#rate-window) hid the drift, but left the rate wherever the first step landed.

## Cause

The cause was not the model or the coalescent approximation. It was the FFT convolution of two functions, `fn convolution_function_function_fine()` in [packages/treetime-distribution/src/distribution_ops/convolve.rs](../../packages/treetime-distribution/src/distribution_ops/convolve.rs). It computed a Riemann sum of the samples, which has a first-order error wherever an operand has a large density at a grid end ([M-distribution-fft-convolution-hard-edge-first-order-error.md](../issues/M-distribution-fft-convolution-hard-edge-first-order-error.md)). In the backward pass, both operands of most edge convolutions have such an end:

- A branch without substitutions has its maximum at zero length
- The coalescent prior factor rises toward the tips (about 2 per year here), so an internal node's subtree distribution peaks at the date of its earliest child

Grid spacings scale with $1/\mu$, so the error per internal node scales with $1/\mu$ as well. With the prior, it was about 0.3 nats per node at $\mu=10^{-4}$ and 0.15 at $3\cdot10^{-4}$. This is one source of the bias in [H-timetree-marginal-likelihood-of-node-times-biased-at-300-grid-points.md](../issues/H-timetree-marginal-likelihood-of-node-times-biased-at-300-grid-points.md). The error is set by the spacing of the branch grids, so more message grid points did not reduce it.

## Evidence

On `const/rep00` (400 tips), every subtree log mass was compared with an independent fine-grid backward pass (Python, two resolutions, extrapolated) that used identical peak-normalized Poisson branch likelihoods and the node factors of TreeTime's prior. Root log mass error, TreeTime minus reference:

| | $\mu=10^{-4}$ | $\mu=3\cdot10^{-4}$ | slope error per unit $\ln\mu$ |
| --- | --- | --- | --- |
| without prior | +9.2 | +9.5 | −3 to −8 |
| prior, Riemann sum | +121 | +61 | −55 |
| prior, 10× finer branch grid | +11.6 | +4.2 | −7 |
| prior, trapezoid rule | −1.7 | +2.2 | about −3, within the reference precision |

An MCMC over the 399 node times with the exact Kingman density and Poisson branch likelihoods puts the maximum of the marginal likelihood at 0.97–1.39 times the true rate in 7 replicates. The mean-field posterior of TreeTime was already close to the exact posterior: its $\sum k - \mu L\,\mathbb E[\sum\tau]$ matched the exact slope (+9.3 against +7.6 at $3\cdot10^{-4}$). Only the computed $\log Z$ was wrong.

Ruled out along the way: the normalization of the mean-field prior (useful for comparing runs, `CoalescentPrior::log_normalization()`, but rate-neutral), the linearization of the prior, the branch-length model, and the number of message grid points.

## Fix

`fn trapezoid_from_overlap_sums()` in [packages/treetime-ops/src/convolution.rs](../../packages/treetime-ops/src/convolution.rs) integrates each overlap by the trapezoid rule, for every convolution algorithm of `treetime-ops`. The result converges at second order (`test_convolution_function_function_hard_edges_converge_at_second_order`). v0 takes the same Riemann sum ([packages/legacy/treetime/treetime/node_interpolator.py#L214-L227](../../packages/legacy/treetime/treetime/node_interpolator.py#L214-L227)), so this is a divergence from v0: [kb/decisions/distribution-convolution-trapezoid-rule.md](../decisions/distribution-convolution-trapezoid-rule.md).

## Remaining problems and options

- **The rate step misses the maximum.** After the fix, the profile of $\log Z$ peaks near the true rate but has a flat top (within a few nats over a factor of 1.5). The three-probe parabola overshoots to the 2-fold bound in the first round, and the later interior window (±25%) then holds the rate away from the peak. Final rates in 7 replicates were 0.52–1.99 times the true rate, against 0.97–1.39 for the exact maximum. Following the slope every round oscillates instead. Options:
  - A least-squares quadratic fit over 5–7 probes spanning about ±60% in $\ln\mu$
  - A bounded 1D search on the profile, now that the profile no longer drifts
  - Revisit the rate window, which was chosen while $\log Z$ was biased
- **Residual few-nat noise in $\log Z$.** Other first-order steps may remain: coarsening of convolution results, mass-window cuts, soft-tail fits. Options:
  - Exact convolution of exponential edge segments ([M-distribution-convolution-drops-operand-tails.md](../issues/M-distribution-convolution-drops-operand-tails.md))
  - A grid spacing near a hard edge tied to the local decay scale
  - Check the other steps for second-order convergence one by one
- **Earlier evidence measured with the bias.** Rate and date results in [timetree-clock-rate-by-marginal-likelihood.md](../decisions/timetree-clock-rate-by-marginal-likelihood.md) and [timetree-frozen-lineage-counts-for-coalescent-prior.md](../decisions/timetree-frozen-lineage-counts-for-coalescent-prior.md) include it. In particular, the low-signal "tilt of about 77 per unit of $\ln$ rate" was mostly this error, not a property of the frozen prior.

Tools: `test_scripts/poisson_timetree_surface.py` (Poisson + Kingman likelihood, landscape, Hessian, MCMC slope of the marginal likelihood).
