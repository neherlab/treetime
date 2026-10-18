# Marginal node-time posteriors are miscalibrated by the steps of the frozen lineage counts

> [!IMPORTANT]
> **Decision required.** The coalescent prior of the time passes needs a different local approximation of the Kingman
> density. v0 has the same steps. Removing the node's own step and smoothing $k(t)$ were tried before without a
> clear result; the exact comparison below now separates the two errors and can score candidates before they reach
> the pipeline.

## Summary

With a coalescent prior, the marginal node-time posteriors that `--time-marginal` reports are skewed against the
simulated truth, and the skew is caused by the form of the prior, not by the time passes or the branch likelihoods.
The prior evaluates the lineage counts $k(t)$ frozen at point estimates of the node times
([timetree-frozen-lineage-counts-for-coalescent-prior.md](../decisions/timetree-frozen-lineage-counts-for-coalescent-prior.md)),
so each node's factor steps at every estimated merger time, including its own.

## Evidence

### Calibration on simulated data

Validation suite, `const` schedule, replicates 0-9, all eight $(L, \mu)$ sets, true topology and root, true clock
rate, true $T_c = 10$ fixed (`--coalescent 10`), `--time-marginal only-final`; 80 runs, 399 internal nodes each.
PIT is the posterior CDF at the true date; calibrated posteriors give PIT uniform on $[0, 1]$.

| | 50% coverage | 90% coverage | 99% coverage | mean PIT | PIT < 0.1 | PIT > 0.9 |
| --- | --- | --- | --- | --- | --- | --- |
| `--coalescent 10` | 0.47 | 0.87 | 0.98 | 0.585 | 0.068 | 0.170 |
| `--coalescent-opt` | 0.47 | 0.87 | 0.98 | 0.565 | 0.077 | 0.158 |
| calibrated | 0.50 | 0.90 | 0.99 | 0.500 | 0.100 | 0.100 |

- The skew persists at the highest signal ($L = 10^4$, $\mu = 10^{-2}$: mean PIT 0.558, 15% above 0.9), where
  `opt`, `skyline` and the fixed true $T_c$ give identical numbers; without a prior the same cell gives mean PIT 0.60,
  so the prior acts there even though its $T_c$ does not.
- Point estimates are unbiased there (date bias $-0.004$ years); the posterior's old tail is too heavy.
- The sign depends on the node's position relative to the sampling window (tips span 0-3.9 years before the latest):
  nodes estimated 0-3 years old have mean PIT 0.60-0.63 (posterior too old), nodes estimated older than 4 years
  0.38-0.44 (posterior too young). Stratify by estimated, not true, age: conditioning on the truth is not a valid
  calibration test.

### Exact posterior on a three-tip tree

`test_marginal_exact.rs` ([`packages/treetime/src/timetree/inference/__tests__/test_marginal_exact.rs`](../../packages/treetime/src/timetree/inference/__tests__/test_marginal_exact.rs))
integrates the posterior of `((A,B)X,C)R` with tips at 2012, 2010 and 2011 and JC69 branches of known substitution
counts. Without a prior the time passes match it to within 0.02 in CDF (zero-count branches and the hard bound at
B included), so message passing and the branch-likelihood builder are exact.

With a constant-$T_c$ prior, an exploratory extension compared three priors against the exact Kingman posterior,
$\log p(x, r) = -[3(2010 - x) + (x - r)]/T_c$ for $r < x < 2010$. The prior of the time passes telescopes to
$\sum_v \ln \frac{\hat k(t_v)(\hat k(t_v) - 1)}{2 T_c} - \sum_{\text{branches}} \int \frac{\hat k(s) - 1}{2 T_c} ds$
with $\hat k$ frozen (and 2 before the earliest merger), and its posterior, integrated numerically, reproduces the
time passes to 0.003-0.08 in CDF. Mean PIT of the exact posterior (0.5 means agreement), $\hat k$ from the exact
posterior medians:

| counts (A, B, X, C) | $T_c$ | frozen $\hat k$: X / R | own step removed: X / R |
| --- | --- | --- | --- |
| 3, 0, 2, 6 | 1 | 0.43 / 0.49 | 0.51 / 0.37 |
| | 10 | 0.38 / 0.47 | 0.50 / 0.48 |
| 1, 0, 0, 2 | 1 | 0.37 / 0.39 | 0.42 / 0.36 |
| | 10 | 0.34 / 0.36 | 0.42 / 0.40 |
| 8, 4, 5, 12 | 1 | 0.56 / 0.54 | 0.54 / 0.34 |
| | 10 | 0.40 / 0.46 | 0.51 / 0.48 |
| 0, 0, 1, 1 | 1 | 0.40 / 0.46 | 0.49 / 0.41 |
| | 10 | 0.37 / 0.44 | 0.48 / 0.47 |

## Mechanism

Two errors, of different scaling:

- **The node's own step.** The merger factor $\hat k(\hat k - 1)/(2 T_c)$ jumps by $\binom{k+1}{2} / \binom{k}{2}$ at
  the node's own estimated time, independent of $T_c$. It dominates when the prior is weak against the data
  ($T_c = 10$ above); removing it brings three of four cases within 0.03 in CDF of the exact posterior.
- **Frozen $\hat k$ in the survival term.** Moving a node changes the lineage count that the other branches see. The
  exact conditional log-prior of a node time $x$ has slope $(k_{-v}(x) - 1)/T_c$, where $k_{-v}$ counts the lineages
  without the node's own merger; a frozen $\hat k$ gives about half of it. This error scales with $1/T_c$ and
  dominates when the prior is strong ($T_c = 1$), where removing the own step over-corrects the root. A node whose
  posterior straddles a neighbour's step (case 1, 0, 0, 2) stays off at both $T_c$.

The full Kingman density couples all node times through $k(t)$; the linear-time passes need a factor local to each
node, so $k(t)$ enters as a mean field and some approximation remains.

## Options

- Remove each node's own step (needs the node's previously committed time in the backward pass).
- Replace the steps of the other nodes by their marginal CDFs from the previous round:
  $\bar k_{-v}(t) = \sum_{\text{tips}} [t < t_i] - \sum_{u \ne v} P(t_u > t)$.
- Score each candidate with the exact three-tip comparison (and a four- or five-tip extension) before running the
  calibration suite.

## Mean-field prior: first implementation

Branch `claude-meanfield-prior` replaces the prior of the time passes by the mean-field prior of the second option
([`mean_field.rs`](../../packages/treetime/src/coalescent/mean_field.rs)): the expectations come from the prior-free
pass, are rebuilt once from the first prior pass and then frozen. On `((A,B)X,(C,D)Y)R` with overlapping sibling clades
the node marginals match the exact Kingman posterior to within 0.03 in CDF (`test_marginal_exact.rs`), where the frozen
lineage counts were off by 0.08-0.38.

Validation suite, keep root, true topology, against the same code with the frozen lineage counts (calibration: 80-240
runs per row, all nodes pooled; accuracy: per-cell means of `keep-{opt,skyline}`):

| | frozen $\hat k$ | mean field |
| --- | --- | --- |
| true $T_c$ fixed, true rate: 90% coverage / mean PIT / PIT > 0.9 | 0.870 / 0.585 / 0.170 | 0.888 / 0.537 / 0.133 |
| `opt`, true rate: 90% coverage / mean PIT / PIT < 0.1 | 0.868 / 0.561 / 0.084 | 0.822 / 0.414 / 0.213 |
| `skyline`, true rate: 90% coverage / mean PIT / PIT < 0.1 | 0.830 / 0.588 / 0.085 | 0.811 / 0.408 / 0.222 |
| estimated rate, $L = 10^4$, $\mu = 3 \cdot 10^{-4}$: mean rate error / mean $\lvert\log\rvert$ rate error | $-12.1\%$ / 0.131 | $+3.2\%$ / 0.049 |
| same: root error / date RMSE, years | 1.99 / 0.215 | 1.43 / 0.184 |
| true rate, $L = 10^3$, $\mu = 3 \cdot 10^{-4}$: root error / date RMSE, years | 6.14 / 0.539 | 8.60 / 0.642 |
| `opt`, $L = 10^3$, $\mu = 3 \cdot 10^{-4}$, estimated rate: $T_c$ log bias | +0.61 | -0.21 |
| run time, all 2880 runs | 16594 s | 17429 s |

- With the true $T_c$ the posteriors in the sampling window are calibrated (mean PIT 0.50-0.52 by estimated age, was
  0.60-0.63); nodes estimated older than 4 years are now too young (mean PIT 0.64-0.69, about 2100 of 32000 nodes).
- The clock-rate bias of the frozen lineage counts at high signal is gone, and accuracy at $L = 10^4$ is equal or
  better in every cell.
- With $T_c$ estimated, the posteriors overshoot to too young, most at low signal (`opt`, $L = 10^3$,
  $\mu = 3 \cdot 10^{-4}$: mean PIT 0.25), and root error and date RMSE at $L = 10^3$, $\mu \le 10^{-3}$ get worse.
  $T_c$ is fitted to the lineage counts at the point times while the prior uses expectations; with the mean field the
  fitted $T_c$ drops by 0.3-0.8 in log at low signal, which strengthens the prior and pulls the nodes younger.
- `data/ebola/{20,100,362}`: no failures and no oscillation; `opt` converges in 6, 3 and 8 rounds (was 7, 12, 12).

## Tc from expected statistics

The skyline fit, and the constant $T_c$ as its one-segment case, now take the expected statistics under the node
posteriors of the round, for independent node times: $I_i = \int_i E[k(k-1)/2]\,dt$ with
$E[k(k-1)] = \bar k^2 - \bar k + V$, $V = \sum_u w_u^2 S_u (1 - S_u)$, and $M_i = \sum_v m_v P(t_v \in i)$
(`LineageExpectation::segment_statistics`), the M-step of an EM iteration for $T_c$. Same validation, keep root:

| | frozen $\hat k$ | mean field | mean field, EM $T_c$ |
| --- | --- | --- | --- |
| `opt`, true rate: 90% coverage / mean PIT / PIT < 0.1 / PIT > 0.9 | 0.868 / 0.561 / 0.084 / 0.152 | 0.822 / 0.414 / 0.213 / 0.075 | 0.886 / 0.481 / 0.122 / 0.095 |
| `skyline`, true rate: same | 0.830 / 0.588 / 0.085 / 0.194 | 0.811 / 0.408 / 0.222 / 0.075 | 0.883 / 0.480 / 0.124 / 0.095 |
| `opt`, estimated rate: same | 0.863 / 0.583 / 0.072 / 0.170 | 0.831 / 0.426 / 0.199 / 0.081 | 0.889 / 0.500 / 0.108 / 0.106 |
| true rate, $L = 10^3$, $\mu = 3 \cdot 10^{-4}$: root error / date RMSE, years | 6.14 / 0.539 | 8.60 / 0.642 | 4.13 / 0.434 |
| true rate, $L = 10^3$, $\mu = 10^{-3}$: same | 2.07 / 0.281 | 2.66 / 0.317 | 1.68 / 0.231 |
| estimated rate, $L = 10^4$, $\mu = 3 \cdot 10^{-4}$: mean rate error / mean $\lvert\log\rvert$ rate error | $-12.1\%$ / 0.131 | $+3.2\%$ / 0.049 | $-5.5\%$ / 0.070 |
| same: root error / date RMSE, years | 1.99 / 0.215 | 1.43 / 0.184 | 1.36 / 0.171 |
| estimated rate, $L = 10^3$, $\mu = 3 \cdot 10^{-4}$: mean $\lvert\log\rvert$ rate error | 0.499 | 0.563 | 0.653 |
| `skyline` $T_c$ log bias, $L = 10^4$ | $-0.18$ to $-0.24$ | $-0.19$ to $-0.32$ | $-0.02$ to $+0.05$ |
| run time, `keep-{opt,skyline}` | 6039 s | 6410 s | 6521 s |

- Root error and date RMSE are equal or lower than with the frozen $\hat k$ in every cell but one ($L = 10^4$,
  $\mu = 10^{-2}$, true rate: root error 0.255 against 0.246 years).
- Calibration is close to nominal except at the lowest signal ($L = 10^3$, $\mu = 3 \cdot 10^{-4}$: mean PIT 0.40);
  by estimated age the mean PIT is 0.47-0.52 throughout.
- The constant $T_c$ (`opt`) stays about 0.3 high in log against the true $T_c(t)$ of the `exp` and `piecewise`
  schedules, a misspecification; the skyline is now nearly unbiased.
- The estimated rate at the lowest signal got worse. In one run without temporal signal (`reroot-opt-est`,
  `piecewise`, rep09, $L = 10^3$, $\mu = 3 \cdot 10^{-4}$, where the other variants stop on a non-positive starting
  rate) the rate fell to $10^{-5}$ (true $3 \cdot 10^{-4}$) and 46 nodes were written without a date; the run
  exited 0. Not yet investigated.

## Related

- [N-coalescent-smoothing-frozen-lineage-counts.md](N-coalescent-smoothing-frozen-lineage-counts.md): the step
  function as a suspected cause of bistable nodes
- [M-timetree-confidence-interval-deficiencies.md](M-timetree-confidence-interval-deficiencies.md): how the reported
  interval is built from these posteriors
- Without a coalescent prior the posteriors are far too old (mean PIT up to 0.97, 90% coverage 0.18 at $L = 10^3$,
  $\mu = 3 \cdot 10^{-4}$, true rate): the flat measure on node times favours old nodes. With an estimated rate the
  intervals omit the rate uncertainty, and coverage falls to 0.67-0.71 for runs whose rate is off by more than 20%.
