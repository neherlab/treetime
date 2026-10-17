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
  gives a parabola, whose maximum, at most 2-fold away, is the target; without a maximum the target is that bound on
  the side of the better probe. The rate moves half-way to the target in log rate, so at most by a factor of $\sqrt 2$
  per round.
- The step does not optimize precisely: the refinement loop repeats it every round, with the coalescent timescale and
  the node times re-estimated in between. A full Brent search (to `1e-3` within a 3-fold window) made the loop
  oscillate: on a flat `log_lh_time` the rate moved a few percent per round, and node times jumped by up to 1.8 years
  between rounds on a 400-tip tree with true rate and topology. The half step damps that coupling, and its three
  backward passes replace the 10-20 of the search.
- The clock model is then fitted with that rate fixed (intercept and outputs as for `--clock-rate`).
- The pre-loop and the initial round keep the regression, which supplies the starting rate and the rooting.

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

The half step against the full search, 10 rounds, `const`, `exp` and `piecewise` schedules, replicates 0-9 (median
date RMSE and root error; "late jump" is the largest node-time change after round 3, median over sets):

| set | rate step | rate bias, mean log | date RMSE | root error | late jump | run time |
| --- | --- | --- | --- | --- | --- | --- |
| $L=10^4$, $\mu=3\cdot10^{-4}$ | search | +0.024 | 0.176 | 1.32 | 0.94 | 150 s |
| | half step | -0.009 | 0.170 | 1.41 | 0.77 | 78 s |
| $L=10^3$, $\mu=10^{-3}$ | search | +0.043 | 0.366 | 3.99 | 0.73 | 122 s |
| | half step | +0.044 | 0.365 | 3.43 | 0.48 | 48 s |
| $L=10^3$, $\mu=3\cdot10^{-4}$ | search | -0.652 | 0.414 | 5.51 | 0.51 | 107 s |
| | half step | -0.723 | 0.406 | 5.41 | 0.41 | 41 s |

On `data/ebola/20` the half step converges in 6 rounds with `--coalescent-opt` (6 with the search) and in 9 with
`--coalescent 2.56` (3 with the search): halving makes the approach to a fixed point geometric.

## Accepted limitations

- The coalescent prior's lineage counts are frozen
  ([timetree-frozen-lineage-counts-for-coalescent-prior.md](timetree-frozen-lineage-counts-for-coalescent-prior.md)),
  which biases the rate at low signal and makes the curvature of `log_lh_time` an overconfident rate uncertainty.
- The fixed-rate clock model carries no regression covariance, so `--covariation` gives no rate uncertainty for these
  runs; the rate-susceptibility analysis then needs `--clock-std-dev`.

## Related

- [distribution-tails-and-arithmetic.md](distribution-tails-and-arithmetic.md#amplitude-re-peaking-and-normalization):
  the amplitude that carries the marginal likelihood through the backward pass
