# Clock rate step: fixed probe window and uncalibrated rate uncertainty

Phase 4 of [timetree-coalescent-prior-and-clock-rate-streamlining.md](../proposals/timetree-coalescent-prior-and-clock-rate-streamlining.md). The step window is implemented: the rate follows the slope of $\log Z$ out of the probe window only in the first round, and afterwards moves only toward a maximum between the probes (`RateWindow` in [packages/treetime/src/timetree/optimization/marginal_clock_rate.rs#L102-L121](../../packages/treetime/src/timetree/optimization/marginal_clock_rate.rs#L102-L121); evidence in [kb/decisions/timetree-clock-rate-by-marginal-likelihood.md](../decisions/timetree-clock-rate-by-marginal-likelihood.md#rate-window)). The step has no external rate prior. Three items remain open.

## Curvature-scaled probe window

The rate step probes $\log Z$ at $\mu/1.25$ and $1.25\mu$, with the target capped at a factor 2 (`PROBE_FACTOR` and `MAX_TARGET_FACTOR`, [packages/treetime/src/timetree/optimization/marginal_clock_rate.rs#L14-L15](../../packages/treetime/src/timetree/optimization/marginal_clock_rate.rs#L14-L15)). Under the mean-field prior the curvature of $\log Z$ can scale the probe half-width: $h = \min(\ln 1.25, c/\sqrt{-\partial^2 \log Z})$. Test whether this also gives a usable rate uncertainty.

## One source for `rate_std`

`rate_std` comes from `--clock-std-dev` when given, else from the Gaussian approximation of the timetree (`fn timetree_gaussian_clock()`, called at [packages/treetime/src/timetree/pipeline.rs#L377](../../packages/treetime/src/timetree/pipeline.rs#L377); `fn determine_rate_std()` at [packages/treetime/src/timetree/confidence.rs#L199-L205](../../packages/treetime/src/timetree/confidence.rs#L199-L205)).

- The Gaussian `rate_std` undercovers where the rate is biased: $\lvert z\rvert < 1.645$ for 0.21 and 0.00 of the runs at $L=10^3$, $\mu \le 10^{-3}$, against 0.77-0.83 with signal (proposal, phase 0b result)
- The bias, not the width, limits date coverage: `--confidence` raised the 90% coverage only from 0.885 to 0.887 overall and from 0.854 to 0.866 for nodes older than 10 years, while the low-signal set stays at 0.83
- Compare the Gaussian curvature with the curvature of $\log Z$ under the mean-field prior on rate $z$-scores, and report the better calibrated one

The low-signal marginal likelihood is not flat: at $L=10^3$, $\mu=3\cdot10^{-4}$ it slopes toward low rates by about 77 per unit of $\ln$ rate. A log-normal rate prior centred at the Gaussian rate therefore either changes nothing (log std 1 or 2) or overrides the likelihood (log std near 0.1, or the Gaussian rate standard deviation). `--clock-std-dev` sets the rate uncertainty of the intervals only.

## Cheaper rate term of `--confidence`

The rate susceptibility of `--confidence` reruns the time passes at $\mu + \sigma$, $\mu - \sigma$ and $\mu$ (`fn compute_rate_susceptibility()`, [packages/treetime/src/timetree/confidence.rs#L24-L73](../../packages/treetime/src/timetree/confidence.rs#L24-L73)) for little gain. The Gaussian passes give $\partial m_v/\partial\mu$, so the rate term of a node's variance, $(\partial m_v/\partial\mu)^2 \operatorname{var}(\mu)$, is available without the reruns. Compare coverage and run time.

## Follow-up

The expected lineage counts are rebuilt once and then frozen. A rebuild every round, from an init without a time pass, gave the best calibration of all variants at three rounds (90% coverage 0.900, 0.846 for nodes older than 10 years) but the worst low-signal rate error (0.99). Compare it again with the frozen expectation once the rate step is settled.

Validation: rounds to convergence; low-signal rate error; rate $z$-scores and date coverage on the calibration subset.
