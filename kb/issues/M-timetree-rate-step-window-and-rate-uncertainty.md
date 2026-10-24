# Clock rate step: fixed probe window and uncalibrated rate uncertainty

Phase 4 of [timetree-coalescent-prior-and-clock-rate-streamlining.md](../proposals/timetree-coalescent-prior-and-clock-rate-streamlining.md). The step window is implemented: the rate follows the slope of $\log Z$ out of the probe window only in the first round, and afterwards moves only toward a maximum between the probes; the probe half-width halves after each step to a maximum between the probes, down to 1/8 of the start (`struct RateWindow` in [packages/treetime/src/timetree/optimization/marginal_clock_rate.rs#L125-L149](../../packages/treetime/src/timetree/optimization/marginal_clock_rate.rs#L125-L149); evidence in [kb/decisions/timetree-clock-rate-by-marginal-likelihood.md](../decisions/timetree-clock-rate-by-marginal-likelihood.md#rate-window) and [#shrinking-probe-window](../decisions/timetree-clock-rate-by-marginal-likelihood.md#shrinking-probe-window)). The step has no external rate prior. Three items remain open.

## Curvature-scaled probe window

The rate step first probes $\log Z$ at $\mu/1.25$ and $1.25\mu$, with the target capped at a factor 2, and halves the probe half-width after each interior maximum (`PROBE_FACTOR`, `MIN_PROBE_SHRINK` and `MAX_TARGET_FACTOR`, [packages/treetime/src/timetree/optimization/marginal_clock_rate.rs#L14-L16](../../packages/treetime/src/timetree/optimization/marginal_clock_rate.rs#L14-L16)). Under the mean-field prior the curvature of $\log Z$ can scale the probe half-width: $h = \min(\ln 1.25, c/\sqrt{-\partial^2 \log Z})$. Test whether this also gives a usable rate uncertainty.

## One source for `rate_std`

`rate_std` comes from `--clock-std-dev` when given, else from the Gaussian approximation of the timetree (`fn timetree_gaussian_clock()`, called at [packages/treetime/src/timetree/pipeline.rs#L364](../../packages/treetime/src/timetree/pipeline.rs#L364); `fn determine_rate_std()` at [packages/treetime/src/timetree/confidence.rs#L208-L214](../../packages/treetime/src/timetree/confidence.rs#L208-L214)).

- The Gaussian `rate_std` undercovers where the rate is biased: $\lvert z\rvert < 1.645$ for 0.21 and 0.00 of the runs at $L=10^3$, $\mu \le 10^{-3}$, against 0.77-0.83 with signal (proposal, phase 0b result)
- The bias, not the width, limits date coverage: `--confidence` raised the 90% coverage only from 0.885 to 0.887 overall and from 0.854 to 0.866 for nodes older than 10 years, while the low-signal set stays at 0.83
- Compare the Gaussian curvature with the curvature of $\log Z$ under the mean-field prior on rate $z$-scores, and report the better calibrated one

> [!WARNING]
> **Needs review.** The tilt of about 77 per unit of $\ln$ rate was measured with the rectangle sum in the convolution of two time functions, which biased $\log Z$ toward slow rates ([kb/decisions/distribution-convolution-trapezoid-rule.md](../decisions/distribution-convolution-trapezoid-rule.md)). The report [kb/reports/timetree-low-signal-rate-bias-from-convolution-edges.md](../reports/timetree-low-signal-rate-bias-from-convolution-edges.md) attributes most of the tilt to that error, not to the frozen prior. Neither the tilt nor that attribution has been measured on the current code: [M-timetree-low-signal-rate-bias-claims-not-reproduced.md](M-timetree-low-signal-rate-bias-claims-not-reproduced.md).

The low-signal marginal likelihood is not flat: at $L=10^3$, $\mu=3\cdot10^{-4}$ it slopes toward low rates by about 77 per unit of $\ln$ rate. A log-normal rate prior centred at the Gaussian rate therefore either changes nothing (log std 1 or 2) or overrides the likelihood (log std near 0.1, or the Gaussian rate standard deviation). `--clock-std-dev` sets the rate uncertainty of the intervals only.

## Cheaper rate term of `--confidence`

The rate susceptibility of `--confidence` reruns the time passes at $\mu + \sigma$, $\mu - \sigma$ and $\mu$ (`fn compute_rate_susceptibility()`, [packages/treetime/src/timetree/confidence.rs#L24-L79](../../packages/treetime/src/timetree/confidence.rs#L24-L79)) for little gain. The Gaussian passes give $\partial m_v/\partial\mu$, so the rate term of a node's variance, $(\partial m_v/\partial\mu)^2 \operatorname{var}(\mu)$, is available without the reruns. Compare coverage and run time.

## Follow-up

The expected lineage counts are rebuilt once and then frozen. A rebuild every round, from an init without a time pass, gave the best calibration of all variants at three rounds (90% coverage 0.900, 0.846 for nodes older than 10 years) but the worst low-signal rate error (0.99). Compare it again with the frozen expectation once the rate step is settled.

Validation: rounds to convergence; low-signal rate error; rate $z$-scores and date coverage on the calibration subset.
