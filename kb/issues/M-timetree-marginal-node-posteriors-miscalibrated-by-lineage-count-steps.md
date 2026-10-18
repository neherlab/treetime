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

## Related

- [N-coalescent-smoothing-frozen-lineage-counts.md](N-coalescent-smoothing-frozen-lineage-counts.md): the step
  function as a suspected cause of bistable nodes
- [M-timetree-confidence-interval-deficiencies.md](M-timetree-confidence-interval-deficiencies.md): how the reported
  interval is built from these posteriors
- Without a coalescent prior the posteriors are far too old (mean PIT up to 0.97, 90% coverage 0.18 at $L = 10^3$,
  $\mu = 3 \cdot 10^{-4}$, true rate): the flat measure on node times favours old nodes. With an estimated rate the
  intervals omit the rate uncertainty, and coverage falls to 0.67-0.71 for runs whose rate is off by more than 20%.
