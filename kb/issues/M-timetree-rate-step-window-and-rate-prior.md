# Clock rate step: fixed probe window and no external rate prior

Phase 4 of [timetree-coalescent-prior-and-clock-rate-streamlining.md](../proposals/timetree-coalescent-prior-and-clock-rate-streamlining.md).

- The rate step probes $\log Z$ at $\mu/1.25$ and $1.25\mu$ with the target capped at a factor 2
  ([`marginal_clock_rate.rs`](../../packages/treetime/src/timetree/optimization/marginal_clock_rate.rs)). Under the
  mean-field prior the curvature of $\log Z$ is no longer overconfident; scale the probe half-width by it,
  $h = \min(\ln 1.25, c/\sqrt{-\partial^2 \log Z})$, and test whether it also gives a usable rate uncertainty.
- Without temporal signal the marginal likelihood cannot fix the rate. A log-normal prior on $\mu$ adds $\log p(\mu)$
  to each probe; it would also replace the hard stop on a non-positive regression rate when given.

- The Gaussian `rate_std` of 0b undercovers where the rate is biased: $\lvert z\rvert < 1.645$ for 0.21 and 0.00 of
  the runs at $L=10^3$, $\mu \le 10^{-3}$, against 0.77–0.83 with signal (proposal, phase 0b result). The bias, not
  the width, limits date coverage: `--confidence` raised the 90% coverage only from 0.885 to 0.887 overall and
  0.854 to 0.866 for nodes older than 10 years, while the low-signal set stays at 0.83.
- One source for `rate_std`: compare the Gaussian curvature with the curvature of $\log Z$ under the mean-field
  prior on rate $z$-scores, and report the better calibrated one.
- The rate susceptibility of `--confidence` reruns the time passes at $\mu \pm \sigma$ and at $\mu$ for little gain.
  The Gaussian passes give $\partial m_v/\partial\mu$, so the rate term of a node's variance,
  $(\partial m_v/\partial\mu)^2 \operatorname{var}(\mu)$, is available without them; compare coverage and run time.

- At low signal the rate drifts with the number of rounds, whatever the start: median $\lvert\ln\rvert$ rate error at
  $L=10^3$, $\mu=3\cdot10^{-4}$ from 0.55 to 0.82 (prior-free seed) and 0.51 to 0.91 (init phase) going from two to
  three rounds (proposal, phase 1 result). Rounds should improve the rate, not degrade it.
- A rebuild of the expected lineage counts every round, from an init without a time pass, gave the best calibration
  of all variants at three rounds (90% coverage 0.900, 0.846 for nodes older than 10 years) but the worst
  low-signal rate (0.99). Once the rate step is stable, compare it again with the frozen expectation.

- Cause of the low-signal drift (probe log Z at $L=10^3$, $\mu=3\cdot10^{-4}$, true rate $3\cdot10^{-4}$): the
  marginal likelihood itself rises toward lower rates, by about +17 at the lower probe and -9 to -20 at the upper one
  round after round, while at $\mu \ge 10^{-3}$ both probes fall off and the rate converges. The probes keep the
  mean-field expectation fixed in calendar time while the tree stretches, so nodes pushed into the past, where
  $\bar k \approx 1$, cost almost nothing. Stretching the expectation with the probe overcorrects: with $T_c$ fixed the
  rate runs away upward (probe differences of 500-900), and with $T_c$ stretched along (which leaves the coalescent
  measure unchanged but for the fixed tip dates) it still does (100-150; $1.8\cdot10^{-2}$ against a true
  $10^{-2}$). The rate step therefore depends on how the prior follows a probe; settle it after the coalescent
  statistics (phase 2) and the root (phase 3). Candidates: profile $T_c$ within each probe, an explicit log-normal
  rate prior, or the Gaussian rate at low signal.

Validation: rounds to convergence; low-signal rate error with and without a prior; rate $z$-scores and date
coverage on the calibration subset.
