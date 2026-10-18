# Proposal: one coherent coalescent prior and clock-rate update in the timetree loop

## Summary

Since `rust` the timetree loop gained a clock rate from the marginal likelihood of the node times and, step by step,
a coalescent prior that is calibrated against the exact Kingman posterior. The fixes were made where each problem
appeared, so the pipeline now carries two representations of the lineage count $k(t)$, an initial round with a
prior-free pass, a prior pass and two rebuilds, and special cases for the root. This proposal replaces them with one
design:

1. **Initialize from the clock.** Once the tree is rooted and a rate is known, the clock regression implies a time for
   every clade. These times seed the expected lineage counts and the first $T_c$, so the first time pass already runs
   under the prior and the prior-free pass disappears.
2. **One representation of $k(t)$.** The expected lineage counts of the mean-field prior, $\bar k(t)$ with its
   variance, serve the prior, the $T_c$ fit, the reported coalescent likelihood, the skyline confidence band and the
   merger rates of polytomy resolution. The point-time lineage counts and `CoalescentModel` go.
3. **No special root.** The mean-field factor is proper before the earliest merger and removes each node's own step,
   so the root needs neither `LINEAGES_BEFORE_ROOT` nor a separate contribution.
4. **Keep the rate update.** The parabolic step on the log marginal likelihood by message scaling stays; its probe
   window and an optional external rate prior are the open points.

## What changed since `rust`

| area | commits | state |
| --- | --- | --- |
| marginal likelihood of the node times | `193628ba8` amplitude, `e4769d595` kept incoming tails, `ef857c2b4` `log_lh_time`, `ceb8dffdc` brute-force test, `82581789c` time-axis scaling | infrastructure, keep |
| clock rate from `log_lh_time` | `cae77a650` Brent search, `ef7d61f6a` half step, `1816aa82d` 2/3 step | keep, tune window |
| prior patches | `50bf38c52` root proper before the earliest merger, `a35dcaff4` root at posterior median, `8137f664d` lineage counts rebuilt once | superseded by the mean field, except possibly the median root |
| mean-field prior | `727d7bd08`, `cc10de999` | keep |
| $T_c$ from expected statistics | `ce2cfbcf9` | keep |
| numerics | `963c6c931` hard-approach law in boundary cells (branch `claude-hard-cell-interp`) | prerequisite |

Evidence (validation suite, keep root, true topology;
[M-timetree-marginal-node-posteriors-miscalibrated-by-lineage-count-steps.md](../issues/M-timetree-marginal-node-posteriors-miscalibrated-by-lineage-count-steps.md)):

| | `rust` | frozen $\hat k$, 2/3 step | mean field + EM $T_c$ |
| --- | --- | --- | --- |
| `opt`, true rate: 90% coverage / mean PIT | | 0.868 / 0.561 | 0.886 / 0.481 |
| est. rate, $L=10^4$, $\mu=3\cdot10^{-4}$: median rate error / root error, years | $-2.5\%$ / 2.40 | $-12.1\%$ / 1.99 | $-5.5\%$ / 1.36 |
| true rate, $L=10^3$, $\mu=3\cdot10^{-4}$: root error / date RMSE, years | 8.77 / 0.67 | 6.14 / 0.54 | 4.13 / 0.43 |
| exact 4-tip Kingman posterior, max CDF difference | | 0.08-0.38 | $\le 0.03$ |

## Current pipeline (branch `claude-meanfield-prior`)

```text
estimate_initial_clock      regression, reroot
run_pre_loop                ML branch lengths, reroot, clock filter, ML branch lengths, reroot
run_initial_round           time passes WITHOUT prior
                            setup_coalescent: point lineage counts, expectation, Tc (EM, from prior-free posteriors)
                            time passes WITH prior
                            rebuild point lineage counts and expectation; freeze both
run_refinement_loop         per round: Tc refit (EM, live posteriors)
                                       CoalescentModel(point counts, Tc) -> polytomy merger rates
                                       CoalescentPrior(frozen expectation, Tc) -> prior
                                       relax clock, polytomy resolution, reconstruction, time passes,
                                       clock rate step on log_lh_time
                            convergence metric log_lh_coal from point counts at point times
refine_final_times          CoalescentPrior for rate susceptibility and the final marginal round
gather_results              coalescent output; Tc band from the skyline Hessian
```

Incoherences:

- **Two $k(t)$.** The prior and the $T_c$ fit use the expectation; the polytomy merger rates, `log_lh_coal`, the
  reported likelihood, the skyline segment layout and the $T_c$ band use point counts at point times. The point counts
  carry the root patch `LINEAGES_BEFORE_ROOT` and the `max(0.5, k - 1)` floor of the per-lineage rate.
- **A prior-free pass that biases the start.** Its posteriors have an old root and stretched deep branches
  ([timetree-frozen-lineage-counts-for-coalescent-prior.md](../decisions/timetree-frozen-lineage-counts-for-coalescent-prior.md));
  the first $T_c$ is fitted to them, and one rebuild is needed to recover. On the four-tip test the mean field built
  from prior-free posteriors was off by up to 0.12 in CDF at $T_c = 1$, after one rebuild by 0.017.
- **Root special cases.** `LINEAGES_BEFORE_ROOT` exists to keep the old root factor proper; the mean-field factor of
  the root has slope $(\bar k_{-R} - 1)/T_c = 1/T_c$ before every other merger and is proper by itself. The median
  root (`a35dcaff4`) answered a root that flipped between two peaks under the old prior; whether the mean-field
  posterior still has that bimodality is untested.
- **Frozen versus live.** The frozen-lineage-count decision separates a frozen prior from a live statistic. With the
  EM fit both come from the same object, the expectation, built from different posteriors (frozen for the prior, the
  round's for $T_c$).

## Proposed design

### 1. Clock-implied clade times

After rooting, the regression gives a rate $\mu$ (or `--clock-rate`). For every node $v$ the tips $i$ below it with
dates $t_i$, weights $w_i$ and divergences $d_{v,i}$ from $v$ imply

$$\hat t_v = \frac{\sum_i w_i \,(t_i - d_{v,i}/\mu)}{\sum_i w_i}, \qquad
\hat\sigma_v^2 = \frac{\sum_i w_i\,(t_i - d_{v,i}/\mu - \hat t_v)^2}{\sum_i w_i} \cdot \frac{1}{n_v^{\text{eff}}} .$$

The sums are the `ClockSet` aggregates the regression already propagates up the tree (`t_sum`, `d_sum`, `norm`, and
their squares), so this is a readout, not a new pass. The estimate is the least-squares root position of each clade at
the global rate; it does not respect parent-child order, and a clade with heterochronous tips can place $\hat t_v$
after its earliest tip.

For the prior only the counts matter, but they must stay valid: $k(t) \ge 1$ requires every clade time before its
earliest tip. Options:

- **A. Points, clamped.** $\hat t_v \leftarrow \min(\hat t_v, \min_i t_i)$; survival $S_v$ a step. Simplest.
- **B. Gaussian survival.** $S_v(t) = 1 - \Phi((t - \hat t_v)/\hat\sigma_v)$, truncated at the earliest tip. Smooths
  $\bar k$ by the clade's own clock scatter, which is the information the regression has. Recommended.

From these, `LineageExpectation` and the first $T_c$ (EM statistics, or point statistics for option A).

### 2. One loop

```text
initialize        rate (regression or given), clade times, expectation E_0, Tc_0 from E_0
round r = 1..     prior   = CoalescentPrior(E, Tc)
                  times   = time passes under prior           -> posteriors P_r
                  Tc      = EM fit on Expectation(P_r)
                  rate    = 2/3 step toward argmax log_lh_time(rate | prior)
                  E       = Expectation(P_1) after round 1, then frozen (as validated)
                  polytomies, relaxed clock, reconstruction as now
final             as now, with the final prior
```

`run_initial_round` loses its two extra passes and becomes round 1 of the loop. The expectation policy keeps the
validated "rebuild once, then freeze"; a damped per-round update ($E \leftarrow (1-\alpha) E + \alpha\,
\text{Expectation}(P_r)$ on $\bar k$ and the survivals) is a separate experiment, since rebuilding every round is what
made the old loop chase its own prior.

### 3. One coalescent statistics object

`LineageExpectation` becomes the only $k(t)$:

| consumer | today | proposed |
| --- | --- | --- |
| prior on node times | expectation (mean field) | unchanged |
| $T_c$ fit (constant, skyline) | expectation (EM) | unchanged |
| skyline segment layout | point counts' first and last breakpoint | expectation grid range |
| $T_c$ band | Hessian of the point-count objective | Hessian of the expected objective (same `solve_log_tc`) |
| reported coalescent log-likelihood, `log_lh_coal` | point counts at point times | expected complete-data log-likelihood $E[\log p] = -\sum_i I_i/T_{c,i} - \sum_i M_i \ln T_{c,i}$ |
| polytomy merger rates | $(k-1)/(2T_c)$ at point counts | $(\bar k - 1)/(2T_c)$ |

What goes: `CoalescentModel`, `compute_lineage_counts`, `collect_tree_events`, `compute_lineage_count_distribution`
with `LINEAGES_BEFORE_ROOT`, `compute_integral_merger_rate`, `CoalescentEdgeData` and the edge-based total likelihood,
and `CoalescentSetup.lineage_counts`. `optimize_skyline` takes the segment statistics $(I_i, M_i)$ and the segment
range and no longer needs the graph or node times.

### 4. Root

The root is an internal node like any other in the mean-field factor. Remove `LINEAGES_BEFORE_ROOT` with the point
counts. Re-test the median-root commitment (`a35dcaff4`) on `data/ebola/20` under the mean-field prior: if the root no
longer flips between two peaks, return to the peak for consistency with the other nodes; otherwise keep the median and
record why.

### 5. Clock rate

Keep the parabolic step on $\log Z(\mu)$ by rescaling the branch messages, three backward passes per round, 2/3 of the
way to the target. Open points:

- **Probe window.** The probes sit at $\mu/1.25$ and $1.25\mu$ and the target is capped at a factor 2. With the
  mean-field prior the curvature of $\log Z$ is no longer overconfident (the frozen-count limitation in
  [timetree-clock-rate-by-marginal-likelihood.md](../decisions/timetree-clock-rate-by-marginal-likelihood.md)), so the
  probe half-width can follow it: $h = \min(\ln 1.25, c / \sqrt{-\partial^2 \log Z})$ from the previous round. The
  same curvature gives the rate uncertainty that `--confidence` now lacks for these runs.
- **External rate prior.** Without temporal signal the marginal likelihood cannot fix the rate (validation:
  $L=10^3$, $\mu=3\cdot10^{-4}$). A log-normal prior on $\mu$ adds $\log p(\mu)$ to each probe and needs no other
  change; it also replaces the hard stop on a non-positive regression rate when given.

## Phases

Each phase is code, tests and kb in separate commits, with the core validation suite, the calibration runs
(`--time-marginal only-final`, PIT and coverage), `data/ebola/{20,100,362}` and the exact posterior tests as gates.
The baseline for each phase is the previous one.

| phase | content | gate |
| --- | --- | --- |
| 0 | merge `claude-hard-cell-interp` into the mean-field branch | tests; dengue/500 completes |
| 1 | clock-implied clade times (option B), one loop, initial round removed | calibration and accuracy equal or better than mean field + EM; one time pass fewer in run time |
| 2 | one statistics object: merger rates, $T_c$ band, reported likelihood, segment layout from the expectation; point counts removed | tests rewritten against analytic expected statistics; outputs documented as a v0 divergence |
| 3 | root: median versus peak under the mean-field prior | `data/ebola/20` root stable over rounds |
| 4 | rate: curvature-scaled probe window; optional log-normal rate prior | rounds to convergence; low-signal rate error with and without a prior |
| 5 | kb: decisions for the mean-field prior, the EM fit, the clock-implied start; supersede the frozen-lineage-count and root-before-merger notes; close [M-timetree-initial-round-applies-coalescent-prior.md](../issues/M-timetree-initial-round-applies-coalescent-prior.md) | |

## Open questions

- Clade-time widths (option A or B), and whether $\hat\sigma_v$ should include the rate uncertainty of the regression.
- Expectation update: rebuild once and freeze, or a damped update every round.
- Polytomy resolution: the sweep's merger rate is the only consumer that samples, not integrates; check that
  $\bar k$ in place of integer $k$ does not change its event model beyond the rate.
- v0 parity: the reported coalescent likelihood, the initial round and the root handling all diverge from v0 and need
  decision files.
- Actionable items here are to be extracted into issues once the plan is accepted.

## Related

- [M-timetree-marginal-node-posteriors-miscalibrated-by-lineage-count-steps.md](../issues/M-timetree-marginal-node-posteriors-miscalibrated-by-lineage-count-steps.md)
- [H-timetree-empty-distribution-division.md](../issues/H-timetree-empty-distribution-division.md)
- [timetree-frozen-lineage-counts-for-coalescent-prior.md](../decisions/timetree-frozen-lineage-counts-for-coalescent-prior.md)
- [timetree-clock-rate-by-marginal-likelihood.md](../decisions/timetree-clock-rate-by-marginal-likelihood.md)
- [timetree-root-time-at-posterior-median.md](../decisions/timetree-root-time-at-posterior-median.md)
- [coalescent-analytic-tc-optimization.md](../decisions/coalescent-analytic-tc-optimization.md)
- [N-coalescent-smoothing-frozen-lineage-counts.md](../issues/N-coalescent-smoothing-frozen-lineage-counts.md)
