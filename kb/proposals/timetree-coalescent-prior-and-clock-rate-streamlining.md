# Proposal: one coherent coalescent prior and clock-rate update in the timetree loop

## Summary

Since `rust` the timetree loop gained a clock rate from the marginal likelihood of the node times and, step by step,
a coalescent prior that is calibrated against the exact Kingman posterior. The fixes were made where each problem
appeared, so the pipeline now carries two representations of the lineage count $k(t)$, an initial round with a
prior-free pass, a prior pass and two rebuilds, and special cases for the root. This proposal replaces them with one
design:

0. **Replace the clock regression.** The weighted-least-squares machinery with its variance model ("covariation"),
   the node-then-local-split root search and its optimizer options give way to an ordinary-least-squares regression on
   two-way sufficient statistics with an exact $O(N)$ root search, and to Gaussian message passing for everything the
   variance model was for. This has long been the intended direction.
1. **Initialize from the regression and Gaussian propagation.** The root-to-tip regression keeps the rate, the root
   and the outliers; a Gaussian approximation of the branch-time likelihood, propagated at that rate, gives every
   internal node a marginal in closed form. The marginals seed the expected lineage counts and the first $T_c$, so the
   first time pass already runs under the prior and the prior-free pass disappears.
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

### 0. Replace the clock regression

**Today** ([`clock/`](../../packages/treetime/src/clock/), about 5100 lines). The two-way sums already exist: the
backward and forward `ClockSet` messages per edge are the subtree and outside sums of `rtt_ols.py`, and the default
variance model (zero on internal edges, unit weight per tip) is OLS. What differs is below.

- `ClockSet` carries weighted sums ($\sum t$, $\sum t^2$, $\sum d$, $\sum d^2$, $\sum td$, norm) and
  `propagate_averages` integrates a branch of variance $v = f\,\ell + o$ (plus a tip offset). This is Gaussian
  message passing on regression sums; with $v = 0$ it is OLS, the default everywhere except `--covariation`
  (`clock` sets $f = 2/L$ and tip slack 3, `timetree` $1/L$ and 10, v0 differs from both).
- The root search scores every node's `ClockSet`, then optimizes the split only on the edges next to the best node,
  by an 11-point grid (or Brent, golden section; `--branch-split-*`). `reroot/` duplicates the search generically for
  `optimize`.
- Outliers entered the fit through `clock_from_child`; fixed in 0a.
- The covariance is the unscaled $H^{-1}$; `determine_rate_std` uses it only with `--covariation`, while the augur
  node data report $\sqrt{(H^{-1})_{00}}$ as `rate_std` regardless.

**Replacement**, following `rtt_ols.py` and `gaussian_timetree.py` of
[neherlab/gaussian-timetree](https://github.com/neherlab/gaussian-timetree):

| part | new | covers |
| --- | --- | --- |
| regression | OLS of root-to-tip divergence on date from per-node sums $(n, \sum T, \sum T^2, \sum D, \sum D^2, \sum TD)$: a subtree message (post-order) and an outside message (pre-order) per node; an excluded tip contributes no sums | slope, intercept, $\chi^2$, $R^2$ at any root; outliers out of the fit by construction |
| root search | combined sums at offset $a$ on every edge; RSS is a convex quadratic in $a$, so each edge's optimum is closed form; minimum over all nodes and edges in $O(N)$ | `LeastSquares` and `ClockFilter` rerooting; `MinDev` as the same search at fixed rate 0; `Oldest` and `Tips` unchanged |
| outliers | residuals about the fitted line at the chosen root, IQR rule | `--clock-filter` |
| rate and root uncertainty | Gaussian message passing at the regression rate ($\eta = \mu(D+M+N)/D$): rate and root covariance from the Hessian, converted to $\mu$ | `rate_std` for `--confidence` and the augur node data, in place of `--covariation` |
| node marginals | the same Gaussian passes | the start of the timetree loop (section 1) |

**Kept:** `divergence.rs`, `date_constraints.rs`, the clock CSV rows and date sources (`rtt.rs`), the fixed-rate
intercept of the in-loop clock update, the topology moves of rerooting (`reroot/placement.rs`).

**Retired:** `ClockVarianceParams`, the variance-carrying `propagate_averages`, `clock/find_best_root/` with its
grid, Brent and golden-section split optimizers, the duplicate search in `reroot/` (the `optimize` min-dev search
becomes the fixed-rate-0 case of the new one), and the flags `--variance-factor`, `--variance-offset`,
`--variance-offset-leaf`, `--branch-split-*`, `--covariation` and `--tip-slack` (decided: retired, not remapped).
The outlier rule ($\lvert\text{residual}\rvert > k \cdot \text{IQD}$) and its default $k = 3$ stay as they are;
both are easy to adapt later, and keeping them limits the change surface.

**Consequences:**

- Outputs: the clock-model JSON changes its statistics (regression $\chi^2$, $R^2$ and the Gaussian covariance
  instead of the WLS Hessian); the augur `rate_std` becomes a real standard error. The app (results, the
  covariation switch, contracts, JSON schemas) follows; that part is the engineer's.
- Tests: the frozen-value dengue/100 clock golden master and every clock and covariation snapshot in `dev/smoke.toml`
  change; the v0 runner golden masters (`test_gm_runner/`) call `estimate_clock_model_with_reroot_policy` at a fixed
  v0 rate and keep the root, so they move to the new API with unchanged expectations.
- kb: resolves or obsoletes M-clock-covariation-variance-diverges-from-v0,
  M-clock-filter-regression-uses-covariation, M-timetree-covariation-tip-slack-default-differs-from-v0,
  N-clock-regression-all-negative-rate, N-reroot-clock-search-duplicates-generic-module,
  N-reroot-split-optimizer-default-diverges-from-v0, N-reroot-leaf-variance-offset-diverges-from-v0; needs decisions
  for replacing covariation and for the exact per-edge root search (both v0 divergences).

### 1. Initialization: regression for the rate and root, Gaussian propagation for the internal nodes

The start is cheap and is refined by the loop, so it should be robust rather than exact.

- **Rate, root and outliers from the root-to-tip regression**, as now (least-squares reroot, residual clock filter).
  The regression is empirically the more robust on real data with misdated or erroneous tips, which the simulated
  validation data do not contain.
- **Internal-node marginals from Gaussian propagation at that rate.** The Gaussian timetree
  ([neherlab/gaussian-timetree](https://github.com/neherlab/gaussian-timetree), `gaussian_timetree.py`,
  `gaussian_timetree_derivation.md`; sections below refer to the derivation) approximates a branch with $d$ mutations
  by elapsed time $\sim\mathcal N(\nu/\eta, \nu/\eta^2)$, $\nu = d + 1$. With messages $(S, R, W)$ of mean
  $S - R/\eta$ and variance $W/\eta^2$, one upward and one outside pass give every node's marginal
  $(m_v, \sigma_v^2) = (S_v - R_v/\eta,\ W_v/\eta^2)$. The $+1$ in $\nu$ is the hazard of the event that ends a
  branch (a split, or sampling at a tip), so the geometry runs at $\eta = \mu\,(D + M + N)/D$ for the substitution
  rate $\mu$, with $D$ substitutions, $M$ internal and $N$ terminal branches (Sections 12, 15.3-15.4); the
  regression's $\mu$ is converted that way.
- **Widths include the rate uncertainty.** $\sigma_v$ is the width at a fixed rate; the regression's rate variance
  propagates to $\sigma_v^2 + (\partial m_v/\partial\mu)^2 \operatorname{var}(\mu)$.

The marginals seed the expectation, $S_v(t) = 1 - \Phi((t - m_v)/\sigma_v)$ truncated at the earliest tip below $v$
so that $\bar k \ge 1$; their means need not respect parent-child order, and the counts do not need it. The first
$T_c$ is the EM fit on this expectation. The two methods form a hierarchy: the exact Gamma/Poisson branch factor at
small $d$ cannot be had in the Gaussian family (Section 15.7), and the grid-based time passes are that refinement.

**Evidence** (240 validation sets, true topology and root; scripts outside the repository, `gaussian_timetree.py` and
`rtt_ols.py` with `tension=True`):

| | OLS | Gaussian $\hat\mu$ ($D{:}M{:}N$) | Gaussian $\eta$ |
| --- | --- | --- | --- |
| rate, $L=10^4$, $\mu=3\cdot10^{-4}$: median error / mean $\lvert\log\rvert$ error | $-8.8\%$ / 0.238 | $+10.2\%$ / 0.127 | $+161\%$ / 0.98 |
| rate, $L=10^3$, $\mu=3\cdot10^{-4}$: same | $+40\%$ / 0.451 | $+26\%$ / 0.260 | $+1593\%$ / 2.87 |
| non-positive rate | 4 of 240 sets | none | none |

| internal nodes without the root: RMSE (bias), years | Gaussian at the OLS rate | clade regression at the OLS rate | time passes without prior, true rate |
| --- | --- | --- | --- |
| $L=10^4$, $\mu=10^{-2}$ | 0.063 (0.00) | 0.112 | 0.033 (-0.01) |
| $L=10^4$, $\mu=3\cdot10^{-4}$ | 0.28 (-0.04) | 0.58 | 0.31 (-0.22) |
| $L=10^3$, $\mu=10^{-3}$ | 0.39 (-0.01) | 0.81 | 0.97 (-0.72) |
| $L=10^3$, $\mu=3\cdot10^{-4}$ | 0.47 (+0.04) | 1.21 | 3.65 (-2.68) |
| lineage curve $\int\lvert\bar k - k\rvert\,dt$, lineage-years, $L=10^3$, $\mu=3\cdot10^{-4}$ | 65 | 103 | 1043 |

- The Gaussian marginals are unbiased; the prior-free time passes are too old at low signal even at the true rate
  (the flat-measure volume effect), which is what the current initial round starts from.
- The Gaussian root is too young at every rate tried (mean bias $+0.6$ to $+7.6$ years at the OLS rate), where the
  regression's root is unbiased: the root time comes from the regression.
- At the true rate the Gaussian $\pm1.645\sigma$ covers 90% of the true times; at the OLS rate 73-87%, hence the
  rate term in the widths.
- The Gaussian $\hat\mu$ never fails and is less noisy than OLS here, but runs 0.4-26% high; on clean simulated
  data this does not decide against the regression's robustness on real data.

### 2. One loop

```text
initialize        regression: rate, root, outliers; Gaussian propagation at that rate: node marginals -> E_0, Tc_0
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
| 0a (done) | OLS regression on two-way sums, exact root search, outliers out of the fit, IQR filter; `clock` and `timetree` on it | identical rate, intercept and root wherever the old search found the optimum (OLS default); differences only where the local split search missed it or outliers leaked in; validation suite rate and root recovery |
| 0b | Gaussian message passing at the regression rate: node marginals, rate and root covariance; replaces `--covariation` for `rate_std` | node RMSE, coverage and lineage-curve distance reproduce section 1; `--confidence` coverage on the validation data |
| 0c | retire the old machinery and flags; outputs, contracts, app (with the engineer) | smoke snapshots re-baselined with the reason recorded; dengue/100 golden master rebuilt |
| 1 | merge `claude-hard-cell-interp`; seed the expectation and $T_c$ from the Gaussian marginals; one loop, initial round removed | tests; dengue/500 completes; calibration and accuracy equal or better than mean field + EM; one time pass fewer |
| 2 | one statistics object: merger rates, $T_c$ band, reported likelihood, segment layout from the expectation; point counts removed | tests against analytic expected statistics; outputs documented as a v0 divergence |
| 3 | root: median versus peak under the mean-field prior | `data/ebola/20` root stable over rounds |
| 4 | rate: curvature-scaled probe window; optional log-normal rate prior | rounds to convergence; low-signal rate error with and without a prior |
| 5 | kb: decisions for the clock replacement, the mean-field prior, the EM fit, the regression-and-Gaussian start; supersede the frozen-lineage-count and root-before-merger notes; close [M-timetree-initial-round-applies-coalescent-prior.md](../issues/M-timetree-initial-round-applies-coalescent-prior.md) | |

## Open questions

- The variance offset $o$ for the node widths (0 here), and whether the root survival should be anchored at the
  regression's root time instead of the Gaussian root marginal.
- Decided: `--covariation` and `--tip-slack` are retired; the outlier rule and its default stay for now.
- Expectation update: rebuild once and freeze, or a damped update every round.
- Polytomy resolution: the sweep's merger rate is the only consumer that samples, not integrates; check that
  $\bar k$ in place of integer $k$ does not change its event model beyond the rate.
- v0 parity: the reported coalescent likelihood, the initial round and the root handling all diverge from v0 and need
  decision files.
- Issues: 0a done (branch `claude-clock-exact-root`),
  [M-clock-gaussian-rate-uncertainty-and-node-marginals.md](../issues/M-clock-gaussian-rate-uncertainty-and-node-marginals.md) (0b),
  [M-clock-retire-covariation-and-split-optimizers.md](../issues/M-clock-retire-covariation-and-split-optimizers.md) (0c),
  [M-timetree-loop-starts-from-gaussian-marginals.md](../issues/M-timetree-loop-starts-from-gaussian-marginals.md) (1),
  [M-coalescent-one-expected-lineage-statistics.md](../issues/M-coalescent-one-expected-lineage-statistics.md) (2),
  [N-timetree-root-median-under-mean-field-prior.md](../issues/N-timetree-root-median-under-mean-field-prior.md) (3),
  [N-timetree-rate-step-window-and-rate-prior.md](../issues/N-timetree-rate-step-window-and-rate-prior.md) (4).

## Related

- [M-timetree-marginal-node-posteriors-miscalibrated-by-lineage-count-steps.md](../issues/M-timetree-marginal-node-posteriors-miscalibrated-by-lineage-count-steps.md)
- [H-timetree-empty-distribution-division.md](../issues/H-timetree-empty-distribution-division.md)
- [timetree-frozen-lineage-counts-for-coalescent-prior.md](../decisions/timetree-frozen-lineage-counts-for-coalescent-prior.md)
- [timetree-clock-rate-by-marginal-likelihood.md](../decisions/timetree-clock-rate-by-marginal-likelihood.md)
- [timetree-root-time-at-posterior-median.md](../decisions/timetree-root-time-at-posterior-median.md)
- [coalescent-analytic-tc-optimization.md](../decisions/coalescent-analytic-tc-optimization.md)
- [N-coalescent-smoothing-frozen-lineage-counts.md](../issues/N-coalescent-smoothing-frozen-lineage-counts.md)
