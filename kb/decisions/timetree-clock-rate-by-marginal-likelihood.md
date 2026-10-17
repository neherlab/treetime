# Under a coalescent prior the clock rate maximizes the marginal likelihood of the node times

With a coalescent prior and no fixed `--clock-rate`, each refinement round takes the clock rate that maximizes the log
marginal likelihood of the node times (`log_lh_time`) under the round's prior, instead of the slope of the root-to-tip
regression. Without a coalescent prior the regression stays.

**Type**: Intentional divergence from v0, which estimates the rate by root-to-tip regression in every mode.

**v1 location**: `fn optimize_marginal_clock_rate` in
[`marginal_clock_rate.rs`](../../packages/treetime/src/timetree/optimization/marginal_clock_rate.rs), called from
`fn update_clock_model` in [`round.rs`](../../packages/treetime/src/timetree/round.rs).

## Method

- The branch-time distributions are a rate-independent likelihood on a branch-length grid divided by the rate, so a
  candidate rate rescales their time axis (`fn Distribution.scale_time()`) and costs one backward pass.
- Brent's method on the square root of the rate, within the round's rate divided and multiplied by 3, to a relative
  tolerance of `1e-3`. The search does not iterate to convergence: the refinement loop repeats it every round, with the
  coalescent timescale and the node times re-estimated in between.
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

The rate is closer to the truth in 19 of 24 sets and the date RMSE lower in 20 of 24. Rates come out about 5% low on
average, most at the lowest signal. On `data/ebola/20` the loop still converges: 6 rounds with `--coalescent-opt` (4
before), 3 with `--coalescent 2.56`.

## Accepted limitations

- The coalescent prior's lineage counts are frozen
  ([timetree-frozen-lineage-counts-for-coalescent-prior.md](timetree-frozen-lineage-counts-for-coalescent-prior.md)),
  which biases the rate at low signal and makes the curvature of `log_lh_time` an overconfident rate uncertainty.
- The fixed-rate clock model carries no regression covariance, so `--covariation` gives no rate uncertainty for these
  runs; the rate-susceptibility analysis then needs `--clock-std-dev`.

## Related

- [distribution-tails-and-arithmetic.md](distribution-tails-and-arithmetic.md#amplitude-re-peaking-and-normalization):
  the amplitude that carries the marginal likelihood through the backward pass
