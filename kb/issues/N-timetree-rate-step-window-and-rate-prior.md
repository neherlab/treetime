# Clock rate step: fixed probe window and no external rate prior

Phase 4 of [timetree-coalescent-prior-and-clock-rate-streamlining.md](../proposals/timetree-coalescent-prior-and-clock-rate-streamlining.md).

- The rate step probes $\log Z$ at $\mu/1.25$ and $1.25\mu$ with the target capped at a factor 2
  ([`marginal_clock_rate.rs`](../../packages/treetime/src/timetree/optimization/marginal_clock_rate.rs)). Under the
  mean-field prior the curvature of $\log Z$ is no longer overconfident; scale the probe half-width by it,
  $h = \min(\ln 1.25, c/\sqrt{-\partial^2 \log Z})$, and test whether it also gives a usable rate uncertainty.
- Without temporal signal the marginal likelihood cannot fix the rate. A log-normal prior on $\mu$ adds $\log p(\mu)$
  to each probe; it would also replace the hard stop on a non-positive regression rate when given.

Validation: rounds to convergence; low-signal rate error with and without a prior.
