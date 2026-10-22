# Marginal node-time posteriors are miscalibrated at low signal

## Summary

With a coalescent prior, the marginal node-time posteriors that `--time-marginal` reports are close to calibrated against simulated truth, except at low signal. Four problems remain:

- At the lowest signal ($L = 10^3$, $\mu = 3 \cdot 10^{-4}$) the posteriors are too young: mean PIT 0.40
- Old nodes are biased: the posteriors of nodes far before the sampling window are off-centre
- With an estimated rate, the rate error at the lowest signal is larger than with the v0 form of the prior
- When the rate collapses at low signal, some internal nodes get an empty posterior and no date

## Current prior

The prior of the time passes is a mean field: `CoalescentPrior` ([packages/treetime/src/coalescent/mean_field.rs#L32-L92](../../packages/treetime/src/coalescent/mean_field.rs#L32-L92)) evaluates the expected lineage counts $\bar k(t)$ of `LineageExpectation` ([packages/treetime/src/coalescent/mean_field.rs#L295-L302](../../packages/treetime/src/coalescent/mean_field.rs#L295-L302)) instead of lineage counts frozen at point estimates of the node times, which v0 uses. The expectation is seeded from the Gaussian timetree, rebuilt once from the posteriors of the first prior pass ([packages/treetime/src/timetree/init.rs#L97-L131](../../packages/treetime/src/timetree/init.rs#L97-L131)), and then frozen for the prior; the loop rebuilds it only when the topology changes ([packages/treetime/src/timetree/refinement_loop.rs#L63-L65](../../packages/treetime/src/timetree/refinement_loop.rs#L63-L65)). The $T_c$ fit reads the live expectation of each round ([kb/decisions/timetree-frozen-lineage-counts-for-coalescent-prior.md](../decisions/timetree-frozen-lineage-counts-for-coalescent-prior.md)).

The skyline fit, and the constant $T_c$ as its one-segment case, take the expected statistics under the node posteriors of the round, for independent node times: $I_i = \int_i E[k(k-1)/2]\,dt$ with $E[k(k-1)] = \bar k^2 - \bar k + V$, $V = \sum_u w_u^2 S_u (1 - S_u)$, and $M_i = \sum_v m_v P(t_v \in i)$ (`LineageExpectation::segment_statistics`, [packages/treetime/src/coalescent/mean_field.rs#L514](../../packages/treetime/src/coalescent/mean_field.rs#L514)). This is the M-step of an EM iteration for $T_c$.

The full Kingman density couples all node times through $k(t)$. The linear-time passes need a factor local to each node, so $k(t)$ enters as a mean field and some approximation remains.

## Evidence

### Exact posteriors on small trees

`test_marginal_exact.rs` ([packages/treetime/src/timetree/inference/__tests__/test_marginal_exact.rs](../../packages/treetime/src/timetree/inference/__tests__/test_marginal_exact.rs)) integrates exact posteriors with JC69 branches of known substitution counts:

- Without a prior, on `((A,B)X,C)R` with tips at 2012, 2010 and 2011, the time passes match the exact posterior to within $2.9 \cdot 10^{-3}$ in CDF (node R, both tip branches with zero counts), zero-count branches and the hard bound at B included ([packages/treetime/src/timetree/inference/__tests__/test_marginal_exact.rs#L36-L70](../../packages/treetime/src/timetree/inference/__tests__/test_marginal_exact.rs#L36-L70)). Without a prior no mean field enters, so this difference comes from the discretization ([M-timetree-time-pass-posteriors-deviate-from-the-exact-posterior.md](M-timetree-time-pass-posteriors-deviate-from-the-exact-posterior.md))
- With the mean-field prior, on `((A,B)X,(C,D)Y)R` with overlapping sibling clades, the node marginals match the exact Kingman posterior to within $1.62 \cdot 10^{-2}$ in CDF at $T_c = 1$ and $T_c = 10$ (node X, overlapping clades, $T_c = 1$; [packages/treetime/src/timetree/inference/__tests__/test_marginal_exact.rs#L77-L112](../../packages/treetime/src/timetree/inference/__tests__/test_marginal_exact.rs#L77-L112)), with per-node limits between $10^{-3}$ and $10^{-1}$. Lineage counts frozen at point times were off by 0.08-0.38 on the same tree

### Calibration on simulated data

Validation suite, keep root, true topology. PIT is the posterior CDF at the true date; calibrated posteriors give PIT uniform on $[0, 1]$, so mean PIT 0.5, and PIT below 0.1 and above 0.9 for 10% of nodes each. Calibration rows pool all nodes of 80-240 runs; accuracy rows are per-cell means of `keep-{opt,skyline}`. The first column is the v0 form of the prior, with the lineage counts frozen at point times.

| | frozen $\hat k$ (v0 form) | mean field, EM $T_c$ |
| --- | --- | --- |
| true $T_c$ fixed, true rate: 90% coverage / mean PIT / PIT > 0.9 | 0.870 / 0.585 / 0.170 | 0.888 / 0.537 / 0.133 |
| `opt`, true rate: 90% coverage / mean PIT / PIT < 0.1 / PIT > 0.9 | 0.868 / 0.561 / 0.084 / 0.152 | 0.886 / 0.481 / 0.122 / 0.095 |
| `skyline`, true rate: same | 0.830 / 0.588 / 0.085 / 0.194 | 0.883 / 0.480 / 0.124 / 0.095 |
| `opt`, estimated rate: same | 0.863 / 0.583 / 0.072 / 0.170 | 0.889 / 0.500 / 0.108 / 0.106 |
| true rate, $L = 10^3$, $\mu = 3 \cdot 10^{-4}$: root error / date RMSE, years | 6.14 / 0.539 | 4.13 / 0.434 |
| true rate, $L = 10^3$, $\mu = 10^{-3}$: same | 2.07 / 0.281 | 1.68 / 0.231 |
| estimated rate, $L = 10^4$, $\mu = 3 \cdot 10^{-4}$: median rate error / mean $\lvert\log\rvert$ rate error | $-12.1\%$ / 0.131 | $-5.5\%$ / 0.070 |
| same: root error / date RMSE, years | 1.99 / 0.215 | 1.36 / 0.171 |
| estimated rate, $L = 10^3$, $\mu = 3 \cdot 10^{-4}$: mean $\lvert\log\rvert$ rate error | 0.499 | 0.653 |
| `skyline` $T_c$ log bias, $L = 10^4$ | $-0.18$ to $-0.24$ | $-0.02$ to $+0.05$ |
| run time, `keep-{opt,skyline}` | 6039 s | 6521 s |

The true-$T_c$ row needs no $T_c$ fit, so it measures the mean field alone.

- Root error and date RMSE are equal or lower than with the frozen $\hat k$ in every cell but one ($L = 10^4$, $\mu = 10^{-2}$, true rate: root error 0.255 against 0.246 years)
- Calibration is close to nominal except at the lowest signal ($L = 10^3$, $\mu = 3 \cdot 10^{-4}$: mean PIT 0.40). By estimated age the mean PIT is 0.47-0.52 throughout
- With the true $T_c$ fixed, the posteriors in the sampling window are calibrated (mean PIT 0.50-0.52 by estimated age), and nodes estimated older than 4 years have mean PIT 0.64-0.69 (about 2100 of 32000 nodes). With the Gaussian seed of the expectation, nodes older than 10 years have 90% coverage 0.794 and date bias +0.77 (240 sets, `--coalescent-opt`, estimated rate)
- The constant $T_c$ (`opt`) stays about 0.3 high in log against the true $T_c(t)$ of the `exp` and `piecewise` schedules, a misspecification; the skyline is nearly unbiased
- The estimated rate at the lowest signal is worse than with the frozen $\hat k$; without temporal signal the marginal likelihood cannot fix the rate
- In one low-signal run the root-to-tip regression after rerooting collapses the rate to $1.67 \cdot 10^{-5}$ (true $3 \cdot 10^{-4}$) with either form of the prior. The marginal-likelihood step then moves it to $1.03 \cdot 10^{-5}$, where 46 internal nodes get an empty posterior and no date
- Stratify by estimated, not true, age: conditioning on the truth is not a valid calibration test

> [!IMPORTANT]
> **Investigation required.** The table predates two changes to the loop: the Gaussian seed of the expectation and the interior-maximum rate window, which reduced the median $\lvert\ln\rvert$ rate error at $L = 10^3$, $\mu = 3 \cdot 10^{-4}$ from 0.505 to 0.281 ([kb/decisions/timetree-clock-rate-by-marginal-likelihood.md](../decisions/timetree-clock-rate-by-marginal-likelihood.md#rate-window)). Rerun the calibration suite (`--time-marginal only-final`) on the current code for the four problems of the summary: low-signal mean PIT, old-node PIT and bias by estimated age, low-signal rate error, and the number of nodes without a posterior. The direction of the old-node bias also needs this check: by the PIT definition above, a mean PIT of 0.64-0.69 puts most posterior mass before the true date, which makes these posteriors too old, while earlier analysis described them as too young.

A partial rerun on the code of 2026-10-20 (12 sets, `keep-opt-est` with `--time-marginal only-final`, mean of 2 replicates) measured the 90% coverage only: 0.783 at $L = 10^3$, $\mu = 3 \cdot 10^{-4}$ (mean interval width 0.353 years), 0.024 below the code before the forward-posterior trim; no cell fell by more than that. Mean PIT, the old-node bias and the nodes without a posterior are not remeasured.

### Reproduction in the slow test tier

`packages/treetime/src/timetree/__tests__/test_sim_timetree_recovery.rs` (`just test-slow-rs`) runs replicates of the standard timetree cell of the simulation tests (`timetree_cell()` in `packages/treetime/src/test_utils/simulation.rs`): 20 tips in 5 batches over 2 years, constant $T_c = 5$, random GTR, $L = 10^3$, realized branch lengths, exact tip dates, keep root, `--time-marginal only-final`, 2 iterations. Release build:

- With the true $T_c$ and the true rate fixed at $\mu = 3 \cdot 10^{-3}$, the 90% intervals are calibrated: coverage 0.906 over 20446 replicates (mean error +0.0059, 99.8% interval $[+0.0041, +0.0076]$, margin $\pm 0.03$; `test_sim_timetree_date_intervals_cover_the_true_dates`)
- `test_sim_timetree_date_intervals_cover_the_true_dates_at_low_signal` (ignored) reproduces the low-signal undercoverage: the same cell at $\mu = 3 \cdot 10^{-4}$ with an estimated rate and `--coalescent-opt`, `--confidence` off, so the intervals leave out the rate uncertainty, as in the measurement above. Over the 2267 replicates of the test, the pooled 90% coverage is 0.720 (a node without an interval counts as not covered); a 200-replicate pilot gave a mean coverage of 0.712. The median interval width is 0.82 years, but 319 of 43073 intervals are wider than 100 years (the widest $5.2 \cdot 10^{12}$ years). Two replicates have internal nodes without a date (18 nodes in total), the fourth problem of the summary; the test fails on the first of them (seed 3001584). The mean $\ln$ rate error is $-0.128$
- A probe of the study's constant-$T_c$ layout (400 tips, 10 per batch at $\tau = 0, 0.1, \dots, 3.9$, $T_c = 10$, $\mu = 3 \cdot 10^{-4}$, $L = 10^3$, the same flags) gave a coverage of 0.73-0.83 per replicate, mean 0.786 over 8 replicates, 1.8 s per run, which agrees with the 0.783 above. The 20-tip cell already reproduces the undercoverage, so the test uses it
- `test_sim_timetree_estimated_clock_rate_recovers_the_true_rate` (ignored): the cell at $\mu = 3 \cdot 10^{-3}$ with an estimated rate and `--coalescent-opt` gives a mean $\ln(\hat\mu / \mu)$ of $-0.087$ over 12920 replicates (standard deviation 0.654, 99.8% interval $[-0.104, -0.069]$, margin $\pm 0.03$). In a 400-replicate pilot the median $\ln$ error was $-0.024$ and the mean relative error $\hat\mu / \mu - 1$ was $-0.014$, while 5% of the runs had a rate below 0.39 times the true rate: the low tail of collapsed rates, not a shift of the bulk, sets the mean log error. The log of a noisy estimator also has a mean of about minus half its relative variance, so part of the deviation comes from the spread of the estimate on 20 tips
- At the true rate, the constant $T_c$ matches the maximum-likelihood $T_c$ of the true node times (the integral of $k(k-1)/2$ over time divided by the number of mergers): mean $\ln$ ratio $-0.0069$ over 3813 replicates, 99.8% interval $[-0.0128, -0.0009]$ (`test_sim_timetree_constant_tc_estimate_recovers_the_tc_of_the_true_node_times`)

## Related

- [M-timetree-confidence-interval-widened-to-point-estimate.md](M-timetree-confidence-interval-widened-to-point-estimate.md): how the reported interval is built from these posteriors
- [M-timetree-rate-step-window-and-rate-uncertainty.md](M-timetree-rate-step-window-and-rate-uncertainty.md): the rate step and the calibration of `rate_std`
- Without a coalescent prior the posteriors are far too old (mean PIT up to 0.97, 90% coverage 0.18 at $L = 10^3$, $\mu = 3 \cdot 10^{-4}$, true rate): the flat measure on node times favours old nodes. With an estimated rate the intervals omit the rate uncertainty, and coverage falls to 0.67-0.71 for runs whose rate is off by more than 20%
