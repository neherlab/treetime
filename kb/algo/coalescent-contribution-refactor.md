# Coalescent contribution model

The coalescent subsystem represents lineage state, time-varying $T_c$, and cumulative merger hazard in decimal calendar years. One immutable `CoalescentModel` supplies node contributions during inference and endpoint-derived edge costs for likelihood and optimization.

## Calendar-coordinate state

Let $t$ be a decimal calendar year and $P$ the most recent event. The model stores

$$
H(t)=\int_t^P \kappa(s)\,ds,
\qquad
\kappa(t)=\frac{\max(0.5,k(t)-1)}{2T_c(t)},
$$

where $k(t)$ is the lineage count. Tree events are aggregated in ascending calendar order. The oldest tail contains one ancestral lineage; applying event deltas toward the present leaves zero lineages after the newest sample. Calendar right-continuity at merger breakpoints exposes the sampled-tree-side lineage count, equivalent to v0's left-sided evaluation in time-before-present coordinates.

Constants are coordinate-independent. A nonconstant `Distribution` used as $T_c(t)$ is evaluated in decimal calendar years. Coordinate enforcement by a dedicated type remains tracked in [kb/issues/N-coalescent-time-scale-coordinate-not-type-enforced.md](../issues/N-coalescent-time-scale-coordinate-not-type-enforced.md).

## Node-grouped objective

For a node with $m$ children and total merger rate $\lambda(t)$, the negative-log contributions are

- leaf: $-H(t)$;
- non-root internal node: $(m-1)(H(t)-\log\lambda(t))$;
- root: $(m-1)(H(t)-\log\lambda(t))+H(t)$.

These signed local terms are meaningful as one telescoped whole-tree objective. They reproduce v0's leaf, internal, and root grouping while preserving actual parent multiplicity at multifurcations.

The backward pass combines every child message first, then evaluates one role-specific contribution on the completed distribution's coordinates. Point support remains a point; range endpoints and function pivots are retained; empty distributions remain empty. Formula inputs are rejected because applying a contribution must not choose a new grid. A leaf contribution affects only the temporary message sent to its parent, leaving the stored date constraint unchanged across repeated inference passes.

## Endpoint-derived objective

Likelihood and $T_c$ optimization use inferred child and parent calendar dates. For an edge from parent $p$ to child $c$, the survival cost is

$$
H(t_p)-H(t_c).
$$

The parent's merger-density term is divided among its outgoing edges so their sum equals the node-grouped objective. Reversed finite endpoints are rejected rather than clamped. This path does not use sequence-derived branch-distribution modes as elapsed coalescent time.

## Shared implementation

- [`packages/treetime/src/coalescent/coalescent.rs`](../../packages/treetime/src/coalescent/coalescent.rs): `CoalescentModel`, node costs, edge cost, and pointwise distribution application.
- [`packages/treetime/src/coalescent/lineage_dynamics.rs`](../../packages/treetime/src/coalescent/lineage_dynamics.rs): calendar lineage state.
- [`packages/treetime/src/coalescent/integration.rs`](../../packages/treetime/src/coalescent/integration.rs): shared scalar rates and checked cumulative hazard.
- [`packages/treetime/src/coalescent/edge_data.rs`](../../packages/treetime/src/coalescent/edge_data.rs): inferred endpoint collection and total objective.
- [`packages/treetime/src/timetree/inference/backward_pass.rs`](../../packages/treetime/src/timetree/inference/backward_pass.rs): child-first local application.

Fixed-$T_c$, optimized-$T_c$, and skyline paths all bind their candidate $T_c$ to the same model primitives. Skyline quadrature, extrapolation, simplex initialization, and validation policy remain independently tracked.
