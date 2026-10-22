# Coalescent lineage counts are frozen for the prior and live for the statistic

$k(t)$ plays two roles in the timetree loop, and only one of them may track the times being
inferred. As the **prior** on node times it is frozen at the tree entering the optimization loop,
after one rebuild in the initial round.
As the **statistic** that $T_c$ is estimated from and the likelihood evaluated against it is read
from the live tree. Deriving both from the current times makes the prior self-referential and the
loop does not converge.

**Type**: New v1 mechanism. v0 has no equivalent because its loop does not iterate far enough for
the feedback to matter.

**v1 location**:
[`LineageExpectation` and `CoalescentPrior`](../../packages/treetime/src/coalescent/mean_field.rs), seeded in
[`init.rs`](../../packages/treetime/src/timetree/init.rs) and used per round in
[`refinement_loop.rs`](../../packages/treetime/src/timetree/refinement_loop.rs).

> [!WARNING]
> **Needs review.** The frozen-prior, live-statistic split still holds, but the code below describes the point-time
> lineage counts and `CoalescentModel`, both removed in phase 2 of
> [timetree-coalescent-prior-and-clock-rate-streamlining.md](../proposals/timetree-coalescent-prior-and-clock-rate-streamlining.md).
> The prior is now the mean-field prior of the expected lineage counts `LineageExpectation`, seeded in the init phase
> from the Gaussian node marginals and rebuilt once from the posteriors of its time pass, and again after a round in
> which polytomy resolution changed the topology (the new mergers enter as points for that round's pass); the statistic ($T_c$ fit,
> reported likelihood, convergence metric) is the expectation of the current posteriors. Phase 5 supersedes this file.
> [Normalization of `log_lh_time`](#normalization-of-log_lh_time) describes the current code.

## The two roles

- **Prior.** $k(t)$ sets the merger rate the coalescent imposes on node times. It reaches
  inference through `CoalescentModel`, which
  [`propagate_distributions_backward`](../../packages/treetime/src/timetree/inference/backward_pass.rs)
  applies as a per-node factor via `root_contribution`, `internal_contribution` and
  `leaf_contribution`.
- **Statistic.** $k(t)$ is what $T_c$ is estimated from, in
  [`optimize_skyline`](../../packages/treetime/src/coalescent/skyline.rs), and what the reported
  coalescent likelihood is evaluated against, in
  `compute_coalescent_total_lh` (removed in phase 2).

Before this change `compute_coalescent_model` derived both from the current times. Each pass moved
the times, which moved $k(t)$, which moved the prior the next pass was inferred under: the loop
chased a receding target. Isolating the blocks on `data/ebola/20` showed the destabilizer is not
$T_c$ estimation -- holding $T_c$ fixed at a hand-picked constant still failed to settle.

## Structure

The three quantities enter differently, so the model is assembled from them rather than read off a
tree:

| quantity | status | where it comes from |
| --- | --- | --- |
| $k(t)$, `lineage_counts` | fixed | `compute_lineage_counts(graph)` in the initial round, rebuilt once, before the loop |
| $T_c$, `tc` | estimated | `optimize_skyline` against the live tree, re-run each round when the mode optimizes it |
| $H(t)=\int\kappa$, `expected_mergers` | compound | `compute_integral_merger_rate(tc, lineage_counts)`, materialised at construction |

`CoalescentModel::new(&lineage_counts, &tc)` is the only constructor. `run_timetree` takes
`Option<&CoalescentModel>` rather than a $T_c$ it turns into a model internally, so the caller
decides which $k(t)$ the prior gets. `None` means the run carries no coalescent prior; the absent
option *is* the gate, rather than a model plus a flag telling the function to ignore it.

## One rebuild, then frozen

The prior-free times of the initial round give $k_0$, which seeds the first pass under the prior.
$k(t)$ is then rebuilt once from the times inferred under that prior and frozen. $k_0$ comes from a
tree without a coalescent prior, typically with a much older root and stretched deep branches, so
the prior it imposes disagrees with the times it produces. Rebuilding every round brings back the
self-referential loop.

On simulated data (10 rounds, 30 sets per cell, true topology and root, `--coalescent-opt`), the
rebuild lowers the median root date error from 3.43 to 1.87 years at $L=10^3$, $\mu=10^{-3}$ and
from 1.41 to 0.79 years at $L=10^4$, $\mu=3\cdot10^{-4}$, and the node date RMSE from 0.365 to
0.261 and from 0.170 to 0.153 years. It lowers the clock rate by about 8%. At $L=10^3$,
$\mu=3\cdot10^{-4}$, where 47 substitutions fall on 400 tips, the deep part of the tree then drifts
between rounds, by 3-70 years in 17 of 28 sets, against below one year with $k_0$; the node date
RMSE rises from 0.406 to 0.473 years. These numbers include committing the root at its posterior
median (see [timetree-root-time-at-posterior-median.md](timetree-root-time-at-posterior-median.md)).
Rebuilding every round made the root swing by decades at that signal and gave no better dates elsewhere.

## Freezing, not re-freezing

After the rebuild $k(t)$ is never updated, including across polytomy resolution, which genuinely
adds internal nodes. A skyline $T_c(t)$ has enough expressiveness to absorb that; a single
constant $T_c$ has less. Re-freezing on topology change is the obvious refinement if the fixed-$T_c$
modes turn out to drift, and is a one-line change in `refine_topology`.

## Before the root

$k(t)$ ends at the root of the tree it was built from, but the root of a later pass may move older. A root
placed before the earliest merger has not merged yet, so the count there is 2, and $H(t)$ continues at their
per-lineage rate $1/(2T_c)$ at the earliest $T_c$. A cherry's root prior is then the Kingman density
$(1/T_c)\,e^{-(t_{tip}-t)/T_c}$ on both sides of the reference root.

v0 counts one lineage there and keeps $H$ flat, so the survival penalty stops at the reference root and the
prior on an older root is flat and improper. Integrated over that flat region, the root adds about
$\ln(1/\mu)$ to `log_lh_time` and pulls the marginal clock rate toward slow clocks. With frozen $k_0$ the
reference root is the deep prior-free one and few runs reach the region; with $k(t)$ rebuilt from a prior
pass the inferred root moved decades into it.

From the root on, the node contributions still match v0 (`test_gm_coalescent_model_matches_v0_node_contributions`
compares that range).

## Polytomy resolution uses the same model

Polytomy resolution samples mergers at the per-branch coalescent rate, so it needs a model even for
a run that asked for no coalescent prior. When the coalescent is enabled it gets the identical
object the prior does -- frozen $k_0$, the round's $T_c$ -- so the rate the sampler draws against is
consistent with the prior the times are inferred under. When the coalescent is disabled, `prior` is
`None` while the model still exists, built from a constant $T_c$ estimated by the same one-segment
analytic solve `--coalescent-opt` uses. That replaces v0's `dummy_coalescent_rate`, which is
calibrated per polytomy to the very time window the sampled history has to fit into; see
[timetree-stochastic-polytomy-resolution.md](timetree-stochastic-polytomy-resolution.md).

## The reported coalescent likelihood is a diagnostic

`log_lh_coal` is evaluated with live $k(t)$ while the times were inferred under frozen $k_0$, so it
is not the objective the loop maximizes and is not guaranteed monotone. Its *stability* across
rounds is the useful signal. This is deliberate: the statistic should describe the tree you
actually have.

## Normalization of `log_lh_time`

The mean-field prior applies one factor per node, built from the frozen expected lineage counts $\bar k(t)$. The product of these factors leaves out two terms of the linearized Kingman density that do not depend on the node times:

$$
P = \int \tfrac12\Big(\bar k(t)^2 - \sum_v w_v^2\, \bar S_v(t)^2\Big)\frac{dt}{T_c(t)}, \qquad R = \sum_v w_v^2\, H\big(t_v^{\text{first}}\big)
$$

where:

- $v$: every node with a time, leaves and internal nodes
- $w_v = 1 - c_v$: the signed weight of node $v$ with $c_v$ children ($+1$ for a leaf, $-1$ for a binary internal node, $1 - m$ for a polytomy with $m$ children)
- $\bar S_v(t)$: the survival function of node $v$ in the frozen expectation (a step function for a leaf at its date)
- $T_c(t)$: the coalescent timescale
- $H(t) = \int_{g_0}^{t} T_c^{-1}\,dt'$, with $g_0$ the first point of the lineage grid
- $t_v^{\text{first}}$: the first knot of the survival function of $v$

`CoalescentPrior::log_normalization()` computes $P + R$ once per prior ([packages/treetime/src/coalescent/mean_field.rs](../../packages/treetime/src/coalescent/mean_field.rs)). `fn log_marginal_likelihood()` adds it to the log mass of the root's subtree distribution when the run has a coalescent prior ([packages/treetime/src/timetree/inference/result.rs#L20-L31](../../packages/treetime/src/timetree/inference/result.rs#L20-L31)); in negative log terms, $P + R$ is subtracted from the sum of the node factors. The reported `log_lh_time` (log line and tracelog column) and the rate probes include it. The constant cancels in the posteriors and in the rate step, so inference does not change.

Reason: with $P + R$, the sum of the node factors at point-mass posteriors equals the Kingman negative log density $M \ln T_c + \int \binom{k(t)}{2} \frac{dt}{T_c}$ for a constant $T_c$, where $M$ is the number of mergers (a polytomy with $m$ children counts as $m - 1$) and $k(t)$ the number of lineages. Without it, `log_lh_time` carries an offset that depends on where the lineage grid starts, so on the rate at which the counts were frozen, and runs that froze their counts at different rates report values that cannot be compared. Derivation: [kb/algo/coalescent-contribution-refactor.md](../algo/coalescent-contribution-refactor.md#normalization).

Tests:

- `test_mean_field_normalized_prior_at_its_point_times_is_the_kingman_density` (`packages/treetime/src/coalescent/__tests__/test_mean_field.rs`): a binary tree and a polytomy
- `test_prop_mean_field_normalized_point_prior_matches_the_kingman_density` (`packages/treetime/src/coalescent/__tests__/test_prop_mean_field_normalization.rs`): simulated trees, within $10^{-9}$. Without the normalization it fails; the smallest failing input is two tips at time 0 with $T_c = 0.5$, where the Kingman value is $-0.4961$ and the factor sum $0.0951$

v0 has no mean-field prior, so no v0 value of `log_lh_time` exists to compare with.

## Results

`data/ebola/20`, per-round maximum node-time change in years:

| round | `--coalescent-opt`, live $k(t)$ | frozen $k_0$ | `--coalescent 2.56`, live $k(t)$ | frozen $k_0$ |
| --- | --- | --- | --- | --- |
| 1 | 0.287 | 0.294 | 0.291 | 0.306 |
| 2 | 0.181 | 0.297 | 0.177 | 0.308 |
| 3 | 0.055 | 0.038 | 0.100 | 0.043 |
| 4 | 0.155 | 0.015 | 0.115 | 0.036 |
| 5 | 0.146 | **0.0022** converged | 0.035 | **0.0019** converged |
| 6-10 | 0.072-0.184, no decay | -- | 0.045-0.152, no decay | -- |

`log_lh_coal` also stops drifting: 26.5 to 27.1 across the run, against 26.2 to 28.6 under live
$k(t)$. Under live $k(t)$ the loop walked downhill on its own reported total, trading sequence
likelihood for coalescent prior with nothing arbitrating.

As a side effect the rate-susceptibility passes now hold one prior fixed across the three perturbed
rates. Previously each rebuilt $k(t)$ from its own perturbed times, confounding the sensitivity
being measured with the prior's reaction to it.

## Related

- [timetree-clock-constrained-profile-propagation.md](timetree-clock-constrained-profile-propagation.md)
  -- the feedback that made this instability observable
- [coalescent-analytic-tc-optimization.md](coalescent-analytic-tc-optimization.md) -- the solve
  behind the estimated $T_c$
