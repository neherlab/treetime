# Coalescent contribution model

The coalescent subsystem represents lineage state, time-varying $T_c$ and the coalescent statistics in decimal calendar years. One object, the expected lineage counts `LineageExpectation`, carries the lineage state for every consumer: the mean-field prior on node times, the $T_c$ fit, the polytomy merger rates, the reported coalescent log-likelihood and the convergence trace.

## Expected lineage counts

With merger indicators $I_u(t) = [t_u > t]$ and node weights $w_u = 1 - n_{\text{children}}(u)$, the lineage count is $k(t) = \sum_u w_u I_u(t)$. Its expectation under the node-time distributions $S_u(t) = P(t_u > t)$ is

$$
\bar k(t) = \sum_u w_u S_u(t),
\qquad
V(t) = \sum_u w_u^2 S_u(t)\,(1 - S_u(t)),
$$

on a uniform grid over the span of the node distributions. A node with a point time has a step survival function; its time enters the grid as a zero-width pair of points that holds the left and the right limit, so the trapezoid rule integrates the step exactly. The survival functions come from node posteriors (`LineageExpectation::from_posterior`) or, in the init phase, from the Gaussian node-time marginals truncated at each node's earliest tip (`LineageExpectation::from_gaussian`). Nodes on bad branches and nodes without a time count as surviving at all times, as the lineage count counts them at the present.

## Mean-field prior

The Kingman log-density couples all node times through $k(t)$. Linearized in the indicators of the other nodes around their expectations, the terms of node $v$ are local:

$$
-\log f_v(x) = \int^x \frac{w_v(w_v-1)/2 + w_v \bar k_{-v}(t)}{T_c(t)}\,dt + m_v \ln T_c(x),
\qquad \bar k_{-v} = \bar k - w_v S_v,
$$

with $m_v$ the mergers at $v$. `CoalescentPrior::node_contribution` evaluates it; the backward pass combines every child message first, then applies the contribution on the completed distribution's coordinates through `distribution_multiply_by_fn` in negative-log space. A leaf contribution affects only the temporary message sent to its parent.

## Coalescent statistics

For segments $i$ of $T_c$, the expected pair-rate integral and merger count are

$$
I_i = \int_i \tfrac12\left(\bar k^2 - \bar k + V\right) dt,
\qquad
M_i = \sum_v m_v\, P(t_v \in i)
$$

(`LineageExpectation::segment_statistics`). The $T_c$ fit maximizes the expected complete-data log-likelihood

$$
\ell(T_c) = -\sum_i \frac{I_i}{T_{c,i}} - \sum_i M_i \ln T_{c,i}
$$

(with a smoothness penalty on $\ln T_c$ between skyline segments), which is also the reported coalescent log-likelihood and the `log_lh_coal` of the convergence trace (`fn expected_log_likelihood`, `CoalescentTimescale::expected_log_likelihood`). The segments span from where $\bar k$ rises half a lineage above its value before every node, near the root, to the latest tip (`LineageExpectation::span`). Polytomy resolution draws mergers at the per-lineage rate $\max(1/2, \bar k - 1)/(2T_c)$ (`LineageExpectation::merger_rate_schedule`); before the grid, $\bar k$ is its value before every node.

Constants are coordinate-independent. A nonconstant $T_c$ is evaluated in decimal calendar years; enforcing the coordinate by a type remains tracked in [kb/issues/N-coalescent-time-scale-coordinate-not-type-enforced.md](../issues/N-coalescent-time-scale-coordinate-not-type-enforced.md).

## Implementation

- [`packages/treetime/src/coalescent/mean_field.rs`](../../packages/treetime/src/coalescent/mean_field.rs): `LineageExpectation` and `CoalescentPrior`.
- [`packages/treetime/src/coalescent/skyline.rs`](../../packages/treetime/src/coalescent/skyline.rs): the $T_c$ fit (`optimize_skyline`), its confidence band from the Hessian of $\ell$, and `expected_log_likelihood`.
- [`packages/treetime/src/timetree/coalescent_timescale.rs`](../../packages/treetime/src/timetree/coalescent_timescale.rs): the fixed, constant and skyline time scales and their report.
- [`packages/treetime-distribution/src/distribution_ops/multiply_by_fn.rs`](../../packages/treetime-distribution/src/distribution_ops/multiply_by_fn.rs): `distribution_multiply_by_fn`, the log-space pointwise application of a per-time contribution to a distribution.
- [`packages/treetime/src/timetree/inference/backward_pass.rs`](../../packages/treetime/src/timetree/inference/backward_pass.rs): child-first application of the prior.

The point-time lineage counts, the edge-based coalescent likelihood and `CoalescentModel` were removed in phase 2 of [timetree-coalescent-prior-and-clock-rate-streamlining.md](../proposals/timetree-coalescent-prior-and-clock-rate-streamlining.md); the reported likelihood is no longer v0's edge-based total.

## Remaining work

- Enforce the calendar coordinate of nonconstant $T_c$ values in the type system: [kb/issues/N-coalescent-time-scale-coordinate-not-type-enforced.md](../issues/N-coalescent-time-scale-coordinate-not-type-enforced.md).
- Decide and enforce the skyline extrapolation policy: [kb/issues/N-coalescent-skyline-extrapolation-policy-undecided.md](../issues/N-coalescent-skyline-extrapolation-policy-undecided.md).
- Validate the internally derived skyline boundaries: [kb/issues/N-coalescent-skyline-grid-validation-incomplete.md](../issues/N-coalescent-skyline-grid-validation-incomplete.md).
