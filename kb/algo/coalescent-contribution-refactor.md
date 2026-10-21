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

### Normalization

The node factors carry the first-order part of the linearization. Their product leaves out two terms that do not depend on the node times:

- **Linearization constant.** Around the expectations, $S_uS_v \approx \bar S_uS_v + S_u\bar S_v - \bar S_u\bar S_v$, so the pair term of the Kingman exponent is $\sum_v w_v\int \bar k_{-v}S_v/T_c - P$ with $P=\int \tfrac12\big(\bar k^2-\sum_v w_v^2\bar S_v^2\big)/T_c\,dt$. The factors carry the first part; each pair interaction enters the factors of both of its nodes.
- **Own-survival offset.** The integrals of $T_c^{-1}$ and $\bar k/T_c$ start at the first point $g_0$ of the lineage grid, but the own-survival integral of node $v$ starts at the first knot of its survival. Starting it at $g_0$ adds $-w_v^2\,H(t^{\text{first}}_v)$ per node, with $H=\int_{g_0}T_c^{-1}$; the factors leave out the sum $R=\sum_v w_v^2 H(t^{\text{first}}_v)$.

Without $P+R$, the log of the factor product carries an offset of about $N_{\text{internal}}\,H(g_0)$. Before the root, every $\bar S_v=1$ and $\bar k=\sum_v w_v=1$, so the integrand of $P$ is $\tfrac12(1-\sum_v w_v^2)$, which is $-N_{\text{internal}}$ on a binary tree. That cancels the pair term $\tfrac12 w_v(w_v-1)H$ of the internal nodes, and with $P+R$ the result no longer depends on where the integrals start. `CoalescentPrior::log_normalization` computes $P+R$ once per prior. `fn log_marginal_likelihood()` adds it, so `log_lh_time` and the rate probes approximate the Kingman log marginal likelihood to first order in the deviation of the node times from the frozen expectation. Inference is unchanged: the constant cancels in the posteriors and in the rate step.

At point-mass posteriors the linearization is exact, and the normalized product equals the Kingman density (`test_mean_field_normalized_prior_at_its_point_times_is_the_kingman_density`, a binary tree and a polytomy).

Before the normalization, `log_lh_time` depended on the rate at which the lineage counts were frozen. An estimated-rate run starts from the regression rate, a run with `--clock-rate` from the given rate, and the two put $g_0$ at different times. On 30 `const` validation sets, the two runs at the same final rate differed by the medians below; the gap that remains is the second-order dependence on the frozen expectation:

| Set | Median gap before | Median gap after | Largest gap after |
| --- | --- | --- | --- |
| $L=10^4$, $\mu=3\cdot10^{-4}$ | 165 | 0.5 | 5.3 |
| $L=10^3$, $\mu=10^{-3}$ | 461 | 0.4 | 8.1 |
| $L=10^3$, $\mu=3\cdot10^{-4}$ | 985 | 31 | 79 |

The gaps are $\lvert\log Z_{\text{est}}-\log Z_{\text{fixed}}\rvert$ in nats, with `--coalescent 10 --keep-root --clock-filter 0 --max-iter 4` on the `realized` input tree.

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
