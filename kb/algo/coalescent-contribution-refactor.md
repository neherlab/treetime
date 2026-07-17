# Coalescent contribution model

The coalescent subsystem represents lineage state, time-varying $T_c$, and the cumulative expected merger count in decimal calendar years. One immutable `CoalescentModel` supplies node contributions during inference and endpoint-derived edge contributions for likelihood and optimization.

## Refactor status

The contribution refactor is complete. One shared `CoalescentModel` replaces per-node formula closures and contribution maps; merger rates are evaluated with scalar calendar-coordinate arithmetic; the backward pass combines all child messages before applying the role-specific contribution on the completed support; and leaf contributions affect only temporary outgoing messages while endpoint-derived edge contributions supply the whole-tree likelihood. The coalescent path therefore does not use formula/function multiplication or perform time-before-present conversion during contribution evaluation.

## Calendar-coordinate state

Let $t$ be a decimal calendar year and $P$ the most recent event. The model stores

$$
H(t)=\int_t^P \kappa(s)\,ds,
\qquad
\kappa(t)=\frac{\max(0.5,k(t)-1)}{2T_c(t)},
$$

where $k(t)$ is the lineage count. In code, $H(t)$ is the model's `expected_mergers` field, the expected number of merger events a branch experiences up to $t$, and $\lambda(t)$ is `total_merger_rate`. Tree events are aggregated in ascending calendar order. The oldest tail contains one ancestral lineage; applying event deltas toward the present leaves zero lineages after the newest sample. Calendar right-continuity at merger breakpoints exposes the sampled-tree-side lineage count, equivalent to v0's left-sided evaluation in time-before-present coordinates.

Constants are coordinate-independent. A nonconstant `Distribution` used as $T_c(t)$ is evaluated in decimal calendar years. Coordinate enforcement by a dedicated type remains tracked in [kb/issues/N-coalescent-time-scale-coordinate-not-type-enforced.md](../issues/N-coalescent-time-scale-coordinate-not-type-enforced.md).

## Node-grouped objective

For a node with $m$ children and total merger rate $\lambda(t)$, the negative-log contributions are

- leaf: $-H(t)$;
- non-root internal node: $(m-1)(H(t)-\log\lambda(t))$;
- root: $(m-1)(H(t)-\log\lambda(t))+H(t)$.

These signed local terms are meaningful as one telescoped whole-tree objective. They reproduce v0's leaf, internal, and root grouping while preserving actual parent multiplicity at multifurcations. In code they are `leaf_contribution`, `internal_contribution`, and `root_contribution`, matching v0's signed `node_contribution`.

The backward pass combines every child message first, then applies one role-specific contribution on the completed distribution's coordinates through `distribution_apply_neg_log_weight` in the distribution layer. That application accumulates in negative-log space, shifting by the combined minimum of contribution minus log-amplitude, so the contribution's own dynamic range does not underflow the peak to zero. The child-message product that feeds it still runs in plain probability space, so tail underflow upstream of the contribution is a separate open concern ([kb/issues/M-timetree-backward-pass-plain-space-underflow.md](../issues/M-timetree-backward-pass-plain-space-underflow.md)), as is the product's uniform-grid resampling ([kb/issues/M-distribution-product-grid-resolution-diverges-from-v0.md](../issues/M-distribution-product-grid-resolution-diverges-from-v0.md)). Point support remains a point; range endpoints and function pivots are retained; empty distributions remain empty. Formula inputs are rejected because applying a contribution must not choose a new grid. A leaf contribution affects only the temporary message sent to its parent, leaving the stored date constraint unchanged across repeated inference passes.

## Endpoint-derived objective

Likelihood and $T_c$ optimization use inferred child and parent calendar dates. For an edge from parent $p$ to child $c$, the survival cost is

$$
H(t_p)-H(t_c).
$$

The parent's merger-density term is divided among its outgoing edges so their sum equals the node-grouped objective. Reversed finite endpoints are rejected rather than clamped. This path does not use sequence-derived branch-distribution modes as elapsed coalescent time.

## Shared implementation

- [`packages/treetime/src/coalescent/coalescent.rs`](../../packages/treetime/src/coalescent/coalescent.rs): `CoalescentModel`, and the node and edge contributions (`leaf`/`internal`/`root`/`edge_contribution`).
- [`packages/treetime-distribution/src/distribution_ops/log_cost.rs`](../../packages/treetime-distribution/src/distribution_ops/log_cost.rs): `distribution_apply_neg_log_weight`, the log-space pointwise application of a per-time contribution to a distribution.
- [`packages/treetime/src/coalescent/lineage_dynamics.rs`](../../packages/treetime/src/coalescent/lineage_dynamics.rs): calendar lineage state.
- [`packages/treetime/src/coalescent/integration.rs`](../../packages/treetime/src/coalescent/integration.rs): shared scalar rates and the checked cumulative merger integral (`expected_mergers`).
- [`packages/treetime/src/coalescent/edge_data.rs`](../../packages/treetime/src/coalescent/edge_data.rs): inferred endpoint collection and total log-likelihood (`coalescent_log_likelihood`).
- [`packages/treetime/src/timetree/inference/backward_pass.rs`](../../packages/treetime/src/timetree/inference/backward_pass.rs): child-first application via the distribution layer.

Fixed-$T_c$, optimized-$T_c$, and skyline paths all bind their candidate $T_c$ to the same model primitives.

## Remaining work

### Core model

- Complete coalescent event state after topology changes: [kb/issues/H-timetree-coalescent-events-incomplete-after-topology-change.md](../issues/H-timetree-coalescent-events-incomplete-after-topology-change.md).
- Enforce the calendar coordinate of nonconstant $T_c$ values in the type system: [kb/issues/N-coalescent-time-scale-coordinate-not-type-enforced.md](../issues/N-coalescent-time-scale-coordinate-not-type-enforced.md).

### Skyline

- Decide and enforce the extrapolation policy: [kb/issues/N-coalescent-skyline-extrapolation-policy-undecided.md](../issues/N-coalescent-skyline-extrapolation-policy-undecided.md).
- Define an accuracy contract for cumulative-hazard quadrature: [kb/issues/N-coalescent-skyline-quadrature-contract-undecided.md](../issues/N-coalescent-skyline-quadrature-contract-undecided.md).
- Decide the scale and direction of simplex initialization: [kb/issues/N-coalescent-skyline-simplex-initialization-undecided.md](../issues/N-coalescent-skyline-simplex-initialization-undecided.md).
- Validate skyline grid shape and domain before endpoint access: [kb/issues/N-coalescent-skyline-grid-validation-incomplete.md](../issues/N-coalescent-skyline-grid-validation-incomplete.md).
