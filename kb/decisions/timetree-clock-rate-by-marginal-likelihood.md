# Under a coalescent prior the clock rate follows the marginal likelihood of the node times

With a coalescent prior and no fixed `--clock-rate`, each refinement round moves the clock rate one damped step toward
the maximum of the log marginal likelihood of the node times (`log_lh_time`) under the round's prior, instead of taking
the slope of the root-to-tip regression. Without a coalescent prior the regression stays.

**Type**: Intentional divergence from v0, which estimates the rate by root-to-tip regression in every mode.

**v1 location**: `fn step_marginal_clock_rate` in
[`marginal_clock_rate.rs`](../../packages/treetime/src/timetree/optimization/marginal_clock_rate.rs), called from
`fn update_clock_model` in [`round.rs`](../../packages/treetime/src/timetree/round.rs).

## Method

- The branch-time distributions are a rate-independent likelihood on a branch-length grid divided by the rate, so a
  candidate rate rescales their time axis (`fn Distribution.scale_time()`) and costs one backward pass.
- One parabolic step in log rate: `log_lh_time` at the round's rate and at that rate divided and multiplied by 1.25
  gives a parabola. In the first round its maximum, at most 2-fold away, is the target; without a maximum the target
  is that bound on the side of the better probe (`RateWindow::FollowSlope`). From the second round on, only a maximum
  between the probes counts, and otherwise the rate stays (`RateWindow::InteriorMaximum`); see
  [Rate window](#rate-window). The rate moves 2/3 of the way to the target in log rate (`STEP_FRACTION`), so at most
  by a factor of $2^{2/3} \approx 1.59$ in the first round and $1.25^{2/3} \approx 1.16$ after it.
- The step does not optimize precisely: the refinement loop repeats it every round, with the coalescent timescale and
  the node times re-estimated in between. A full Brent search (to `1e-3` within a 3-fold window) made the loop
  oscillate: on a flat `log_lh_time` the rate moved a few percent per round, and node times jumped by up to 1.8 years
  between rounds on a 400-tip tree with true rate and topology. The partial step damps that coupling, and its three
  backward passes replace the 10-20 of the search. See [Step fraction](#step-fraction) for the choice of 2/3.
- The clock model is then fitted with that rate fixed (intercept and outputs as for `--clock-rate`).
- The pre-loop and the initial round keep the regression, which supplies the starting rate and the rooting.
- After the last round, a final time pass infers the node times at the rate the clock model reports, in every
  `--time-marginal` mode, so the output dates and the reported rate agree (`fn refine_final_times()`
  [packages/treetime/src/timetree/pipeline.rs#L372](../../packages/treetime/src/timetree/pipeline.rs#L372)).

## Why only with a coalescent prior

Without a prior the node times are integrated against a flat measure, and each internal node contributes about
$\ln(1/\mu)$: the marginal likelihood always favours slower clocks. On 72 simulated sets 40 runs ended at the lower
search bound. A coalescent prior is a proper prior on the node times and removes this.

## Evidence

Validation suite, `const` and `exp` schedules, replicates 0-2, $L \in \{10^3, 10^4\}$, $\mu \in \{10^{-3}, 3\cdot10^{-3}\}$,
true topology and root, `--coalescent-opt`, `--clock-filter 0`, default number of rounds, 24 sets:

| | rate error, mean abs. log | node date RMSE, years | root date error, years | run time |
| --- | --- | --- | --- | --- |
| root-to-tip regression | 0.145 | 0.227 | 2.49 | 6.5 s |
| marginal likelihood | 0.077 | 0.186 | 1.88 | 9.4 s |

These numbers are for the full Brent search. The rate is closer to the truth in 19 of 24 sets and the date RMSE lower
in 20 of 24. Rates come out about 5% low on average, most at the lowest signal.

The half step (the fraction before 2/3) against the full search, 10 rounds, `const`, `exp` and `piecewise` schedules,
replicates 0-9 (median date RMSE and root error; "late jump" is the largest node-time change after round 3, median over
sets):

| set | rate step | rate bias, mean log | date RMSE | root error | late jump | run time |
| --- | --- | --- | --- | --- | --- | --- |
| $L=10^4$, $\mu=3\cdot10^{-4}$ | search | +0.024 | 0.176 | 1.32 | 0.94 | 150 s |
| | half step | -0.009 | 0.170 | 1.41 | 0.77 | 78 s |
| $L=10^3$, $\mu=10^{-3}$ | search | +0.043 | 0.366 | 3.99 | 0.73 | 122 s |
| | half step | +0.044 | 0.365 | 3.43 | 0.48 | 48 s |
| $L=10^3$, $\mu=3\cdot10^{-4}$ | search | -0.652 | 0.414 | 5.51 | 0.51 | 107 s |
| | half step | -0.723 | 0.406 | 5.41 | 0.41 | 41 s |

On `data/ebola/20` the half step converges in 6 rounds with `--coalescent-opt` (6 with the search) and in 9 with
`--coalescent 2.56` (3 with the search): the partial step makes the approach to a fixed point geometric.

## Step fraction

With step fraction $a$, the log rate moves $x' = x + a\,(T - x)$ toward the target $T$. Between rounds $T$ follows the
rate through the re-estimated times and Tc; near the fixed point $T = x^* + s\,(x - x^*)$, and the distance to $x^*$
shrinks by $1 - a(1 - s)$ per round. A full step ($a = 1$) oscillates where the target overshoots ($s \approx -1$),
which is what the Brent search did at low signal; a smaller step damps that but converges slowly where the target lags
($0 < s < 1$).

Measured slopes $s$ (within-run secant of $\ln T$ on $\ln$ rate, rounds 2-12, step 1/2, `keep-{opt,skyline}-est`, all
schedules, replicates 0-9): median $+0.30$ at $L=10^4$, $\mu=3\cdot10^{-4}$; $+0.11$ at $L=10^4$, $\mu=3\cdot10^{-3}$;
$-0.10$ at $L=10^3$, $\mu=10^{-3}$, with a quarter of the runs below $-0.5$. So the target lags at high signal and
overshoots at low signal, and no single fraction is best for both.

Same runs, 12 rounds (rounds to reach 5% / 2% of the step-1/2 fixed point, median; late jitter is the within-run SD of
$\ln$ rate over rounds 6-12, median):

| set | step | rounds to 5% / 2% | late jitter | late max node-time change, years |
| --- | --- | --- | --- | --- |
| $L=10^4$, $\mu=3\cdot10^{-4}$ | 1/2 | 4 / 6 | 0.008 | 0.26 |
| | 2/3 | 3 / 4 | 0.007 | 0.28 |
| | 3/4 | | 0.007 | 0.26 |
| | 1 | | 0.011 | 0.30 |
| $L=10^3$, $\mu=10^{-3}$ | 1/2 | 6 / 6.5 | 0.118 | 0.64 |
| | 2/3 | 4 / 5.5 | 0.189 | 0.69 |
| | 3/4 | | 0.243 | 0.81 |
| | 1 | | 0.371 | 1.85 |

2/3 reaches the high-signal fixed point about one round sooner than 1/2, at 1.6 times the late rate jitter at low
signal; 3/4 and above bring the low-signal oscillation back. An adaptive fraction $a = 1/(1-s)$ from the secant of the
last two rounds, bounded to $[0.25, 1]$, was about half a round faster than 2/3 at high signal, noisier at low signal
(mean $|\ln$ rate error$|$ after round 2: 0.32 against 0.26), and no more stable on `data/ebola/{20,100,362}`, so it
was not kept. The rate error against the truth after round 2 differs between the three rules by at most 0.02 in
$\ln$ rate at high signal, against 0.11 left after 12 rounds, most of it the lineage-count bias below.

## Rate window

> [!WARNING]
> **Needs review.** The tilt of about 77 per unit of $\ln$ rate was measured with the rectangle sum in the convolution of two time functions, which biased $\log Z$ toward slow rates ([distribution-convolution-trapezoid-rule.md](distribution-convolution-trapezoid-rule.md)). The report [kb/reports/timetree-low-signal-rate-bias-from-convolution-edges.md](../reports/timetree-low-signal-rate-bias-from-convolution-edges.md) attributes most of the tilt to that error, not to the frozen prior. Neither the tilt nor that attribution has been measured on the current code: [kb/issues/M-timetree-low-signal-rate-bias-claims-not-reproduced.md](../issues/M-timetree-low-signal-rate-bias-claims-not-reproduced.md).

At low signal the rate drifted with the number of rounds: at $L=10^3$, $\mu=3\cdot10^{-4}$ the median
$\lvert\ln\rvert$ rate error grew from 0.51 after two rounds to 1.26 after six, and the node date RMSE from 0.49 to
1.1 years. The probes see the mean-field prior frozen in calendar time while the tree stretches, which tilts
`log_lh_time` toward low rates by about 77 per unit of $\ln$ rate, round after round, without an interior maximum
([M-timetree-rate-step-window-and-rate-uncertainty.md](../issues/M-timetree-rate-step-window-and-rate-uncertainty.md)). Following
that slope to the 2-fold bound every round accumulates the tilt.

Moving the rate only toward a maximum between the probes stops the drift, but alone it also stops the first round from
correcting a poor start: on `const/rep03` at $L=10^3$, $\mu=10^{-2}$ the regression starts at half the true rate, the
probes rise monotonically toward it, and the rate stayed at $4.85\cdot10^{-3}$ against a true $10^{-2}$, which tripled
the date RMSE of the high-signal sets. The first round therefore keeps the slope-following window, and the later
rounds use the interior one. A Gaussian-timetree start rate instead of the regression's fixed `rep03` but raised the
date RMSE at $L=10^3$, $\mu=3\cdot10^{-3}$ from 0.155 to 0.191, so the regression start stays.

Gate: 240 validation sets (`const`, `exp`, `piecewise`, replicates 0-9), true topology and root, `--coalescent-opt`,
`--clock-filter 0`; median $\lvert\ln\rvert$ rate error / date RMSE in years / 90% coverage, two and six rounds:

| set | slope every round, 2 rounds | 6 rounds | slope in round 1 only, 2 rounds | 6 rounds |
| --- | --- | --- | --- | --- |
| $L=10^3$, $\mu=3\cdot10^{-4}$ | 0.505 / 0.494 / 0.830 | 1.259 / 1.109 / 0.856 | 0.281 / 0.503 / 0.820 | 0.290 / 0.536 / 0.776 |
| $L=10^3$, $\mu=10^{-3}$ | 0.143 / 0.273 / 0.875 | 0.191 / 0.345 / 0.837 | 0.191 / 0.272 / 0.873 | 0.208 / 0.343 / 0.833 |
| $L=10^3$, $\mu=3\cdot10^{-3}$ | 0.052 / 0.155 / 0.895 | 0.048 / 0.172 / 0.885 | 0.052 / 0.159 / 0.894 | 0.048 / 0.160 / 0.883 |
| $L=10^3$, $\mu=10^{-2}$ | 0.043 / 0.125 / 0.899 | 0.038 / 0.102 / 0.896 | 0.043 / 0.125 / 0.899 | 0.038 / 0.102 / 0.896 |
| $L=10^4$, $\mu=3\cdot10^{-4}$ | 0.058 / 0.160 / 0.897 | 0.060 / 0.219 / 0.867 | 0.062 / 0.158 / 0.898 | 0.060 / 0.173 / 0.872 |
| $L=10^4$, $\mu=10^{-3}$ | 0.022 / 0.098 / 0.896 | 0.025 / 0.093 / 0.897 | 0.022 / 0.098 / 0.896 | 0.025 / 0.093 / 0.897 |
| $L=10^4$, $\mu=3\cdot10^{-3}$ | 0.018 / 0.056 / 0.903 | 0.018 / 0.056 / 0.902 | 0.018 / 0.056 / 0.903 | 0.018 / 0.056 / 0.902 |
| $L=10^4$, $\mu=10^{-2}$ | 0.010 / 0.039 / 0.907 | 0.012 / 0.034 / 0.907 | 0.010 / 0.039 / 0.907 | 0.012 / 0.034 / 0.907 |

Six rounds no longer degrade the low-signal sets, and the high-signal sets ($\mu \ge 3\cdot10^{-3}$) are unchanged.
Costs: the rate error at $L=10^3$, $\mu=10^{-3}$ (0.143 to 0.191 after two rounds, at the same date RMSE), and the
90% coverage at $L=10^3$, $\mu=3\cdot10^{-4}$ after six rounds (0.856 to 0.776), where the intervals of the drifted
runs had been wide enough to cover dates more than a year off.

Full validation suite (2880 runs, default two rounds, no `--time-marginal`): no failures and no run time change;
median rate error at $L=10^3$, $\mu=3\cdot10^{-4}$ 0.505 to 0.281 (`keep-opt-est`), 0.780 to 0.530 (`reroot-opt-est`),
0.783 to 0.453 (`reroot-skyline-est`), at $L=10^3$, $\mu=10^{-3}$ 0.143 to 0.191 (`keep-opt-est`) and 0.189 to 0.166
(`reroot-skyline-est`), elsewhere within 0.01. No date changed: in these runs the output dates came from the rate
before the last step. The final time pass (see [Method](#method)) has since made the dates follow the last step.

## Accepted limitations

- The coalescent prior holds the mean-field expectation of the lineage counts fixed in calendar time between
  rebuilds ([timetree-frozen-lineage-counts-for-coalescent-prior.md](timetree-frozen-lineage-counts-for-coalescent-prior.md)).
  The expectation is seeded from the Gaussian timetree (or from a time pass without a prior when the Gaussian timetree
  is not available), rebuilt once from the posteriors of the first time pass under the prior, and rebuilt again only
  after a round in which polytomy resolution changed the topology
  ([packages/treetime/src/timetree/init.rs#L70-L93](../../packages/treetime/src/timetree/init.rs#L70-L93),
  [packages/treetime/src/timetree/refinement_loop.rs#L71-L77](../../packages/treetime/src/timetree/refinement_loop.rs#L71-L77)).
  At low signal the tree stretches under this fixed prior, which tilts `log_lh_time` toward low rates (see
  [Rate window](#rate-window)) and biases the rate low.
- The step gives no rate uncertainty. `rate_std` comes from `--clock-std-dev` or from the Gaussian approximation of the
  timetree, and it undercovers where the rate is biased
  ([M-timetree-rate-step-window-and-rate-uncertainty.md](../issues/M-timetree-rate-step-window-and-rate-uncertainty.md)).

## Related

- [distribution-tails-and-arithmetic.md](distribution-tails-and-arithmetic.md#amplitude-re-peaking-and-normalization):
  the amplitude that carries the marginal likelihood through the backward pass
