# Coalescent lineage counts are frozen for the prior and live for the statistic

$k(t)$ plays two roles in the timetree loop, and only one of them may track the times being
inferred. As the **prior** on node times it is frozen at the tree entering the optimization loop.
As the **statistic** that $T_c$ is estimated from and the likelihood evaluated against it is read
from the live tree. Deriving both from the current times makes the prior self-referential and the
loop does not converge.

**Type**: New v1 mechanism. v0 has no equivalent because its loop does not iterate far enough for
the feedback to matter.

**v1 location**:
[`compute_lineage_counts`](../../packages/treetime/src/coalescent/lineage_counts.rs),
[`CoalescentModel`](../../packages/treetime/src/coalescent/coalescent.rs), assembled per round in
[`pipeline.rs`](../../packages/treetime/src/timetree/pipeline.rs).

## The two roles

- **Prior.** $k(t)$ sets the merger rate the coalescent imposes on node times. It reaches
  inference through `CoalescentModel`, which
  [`propagate_distributions_backward`](../../packages/treetime/src/timetree/inference/backward_pass.rs)
  applies as a per-node factor via `root_contribution`, `internal_contribution` and
  `leaf_contribution`.
- **Statistic.** $k(t)$ is what $T_c$ is estimated from, in
  [`optimize_skyline`](../../packages/treetime/src/coalescent/skyline.rs), and what the reported
  coalescent likelihood is evaluated against, in
  [`compute_coalescent_total_lh`](../../packages/treetime/src/coalescent/total_lh.rs).

Before this change `compute_coalescent_model` derived both from the current times. Each pass moved
the times, which moved $k(t)$, which moved the prior the next pass was inferred under: the loop
chased a receding target. Isolating the blocks on `data/ebola/20` showed the destabilizer is not
$T_c$ estimation — holding $T_c$ fixed at a hand-picked constant still failed to settle.

## Structure

The three quantities enter differently, so the model is assembled from them rather than read off a
tree:

| quantity | status | where it comes from |
| --- | --- | --- |
| $k(t)$, `lineage_counts` | fixed | `compute_lineage_counts(graph)`, once, before the loop |
| $T_c$, `tc` | estimated | `optimize_skyline` against the live tree, re-run each round when the mode optimizes it |
| $H(t)=\int\kappa$, `expected_mergers` | compound | `compute_integral_merger_rate(tc, lineage_counts)`, materialised at construction |

`CoalescentModel::new(&lineage_counts, &tc)` is the only constructor. `run_timetree` takes
`Option<&CoalescentModel>` rather than a $T_c$ it turns into a model internally, so the caller
decides which $k(t)$ the prior gets. `None` means the run carries no coalescent prior; the absent
option *is* the gate, rather than a model plus a flag telling the function to ignore it.

## Freezing, not re-freezing

$k(t)$ is frozen once and never updated, including across polytomy resolution, which genuinely
adds internal nodes. A skyline $T_c(t)$ has enough expressiveness to absorb that; a single
constant $T_c$ has less. Re-freezing on topology change is the obvious refinement if the fixed-$T_c$
modes turn out to drift, and is a one-line change in `refine_topology`.

## Polytomy resolution uses the same model

Polytomy resolution samples mergers at the per-branch coalescent rate, so it needs a model even for
a run that asked for no coalescent prior. When the coalescent is enabled it gets the identical
object the prior does — frozen $k_0$, the round's $T_c$ — so the rate the sampler draws against is
consistent with the prior the times are inferred under. When the coalescent is disabled, `prior` is
`None` while the model still exists, built from a constant $T_c$ estimated by the same one-segment
analytic solve `--coalescent-opt` uses. That replaces v0's `dummy_coalescent_rate`, which is
calibrated per polytomy to the very time window the sampled history has to fit into; see
[timetree-stochastic-polytomy-resolution.md](../proposals/timetree-stochastic-polytomy-resolution.md).

## The reported coalescent likelihood is a diagnostic

`log_lh_coal` is evaluated with live $k(t)$ while the times were inferred under frozen $k_0$, so it
is not the objective the loop maximizes and is not guaranteed monotone. Its *stability* across
rounds is the useful signal. This is deliberate: the statistic should describe the tree you
actually have.

## Results

`data/ebola/20`, per-round maximum node-time change in years:

| round | `--coalescent-opt`, live $k(t)$ | frozen $k_0$ | `--coalescent 2.56`, live $k(t)$ | frozen $k_0$ |
| --- | --- | --- | --- | --- |
| 1 | 0.287 | 0.294 | 0.291 | 0.306 |
| 2 | 0.181 | 0.297 | 0.177 | 0.308 |
| 3 | 0.055 | 0.038 | 0.100 | 0.043 |
| 4 | 0.155 | 0.015 | 0.115 | 0.036 |
| 5 | 0.146 | **0.0022** converged | 0.035 | **0.0019** converged |
| 6-10 | 0.072-0.184, no decay | — | 0.045-0.152, no decay | — |

`log_lh_coal` also stops drifting: 26.5 to 27.1 across the run, against 26.2 to 28.6 under live
$k(t)$. Under live $k(t)$ the loop walked downhill on its own reported total, trading sequence
likelihood for coalescent prior with nothing arbitrating.

As a side effect the rate-susceptibility passes now hold one prior fixed across the three perturbed
rates. Previously each rebuilt $k(t)$ from its own perturbed times, confounding the sensitivity
being measured with the prior's reaction to it.

## Related

- [timetree-clock-constrained-profile-propagation.md](timetree-clock-constrained-profile-propagation.md)
  — the feedback that made this instability observable
- [coalescent-analytic-tc-optimization.md](coalescent-analytic-tc-optimization.md) — the solve
  behind the estimated $T_c$
- [N-coalescent-smoothing-frozen-lineage-counts.md](../issues/N-coalescent-smoothing-frozen-lineage-counts.md)
  — $k_0$ is a step function and is not smoothed
