# Implement DivStats divergence-only scoring

Implement `DivStats`, a `RootStats` implementation that scores root positions by variance of root-to-tip distances. No date information required. Implement the two-pass traversal that computes per-edge `(to_parent, to_child)` DivStats from branch lengths.

Depends on [reroot-create-module-with-root-stats-trait.md](reroot-create-module-with-root-stats-trait.md).

## Scope

### DivStats struct (`reroot/div_stats.rs`)

```rust
pub struct DivStats {
    count: f64,
    d_sum: f64,
    dsq_sum: f64,
}
```

`RootStats` implementation:

- `leaf`: `count = 1/var`, `d_sum = bl/var`, `dsq_sum = bl^2/var`. Always contributes (no date gate -- this is the key difference from `ClockSet.norm`).
- `propagate`: `ClockSet::propagate_averages` formulas with all time terms zeroed out. See derivation in the proposal.
- `score`: $(S_{dd} \cdot n - S_d^2) / (2n)$ where $n$ = `count`. Weighted variance of root-to-tip distances. `count` is always positive, so no division by zero.
- `Add`: elementwise sum of all three fields.

### Two-pass traversal (`reroot/div_stats_traversal.rs`)

Compute `BTreeMap<GraphEdgeKey, (DivStats, DivStats)>` from a graph with branch lengths:

1. Leaves-to-root pass: for each edge, accumulate child-side statistics into `to_parent`
2. Root-to-leaves pass: for each edge, compute `to_child` from parent's aggregated statistics minus the child subtree's contribution

Trait bounds: `N: GraphNode, E: GraphEdge` (plus `HasBranchLength` or equivalent). No `ClockNode`/`ClockEdge` required.

Variance model: `branch_variance = variance_factor * bl + variance_offset`, matching `ClockParams` defaults (`variance_factor = 1.0`, `variance_offset = 0.0`). Parameterize via a config struct or use `ClockParams` directly.

## Locations

- `packages/treetime/src/reroot/div_stats.rs`
- `packages/treetime/src/reroot/div_stats_traversal.rs`

## Tests

### Unit tests

- Star tree (n tips, equal branch lengths d): `DivStats.score() = 0` at the true center, `> 0` at any other position
- Two-tip tree: minimum at x = 0.5 (midpoint). Score is a quadratic in x with analytical minimum
- Asymmetric tree (one long branch, one short): minimum biased toward the subtree with more tips
- `DivStats::propagate` cross-check: verify against `ClockSet::propagate_averages` with `t_sum = 0, tsq_sum = 0, dt_sum = 0` on the same inputs

### Integration tests

- `find_best_root::<DivStats>` on flu/h3n2/20: finds a root, score is finite, edge key and split are valid
- Synthetic tree with known optimal root (balanced tree): verify search finds the analytically correct position

### Property tests

- For any tree, `DivStats.score()` at the optimal root <= score at any other position (verified by comparing optimal score against scores at 10 random positions)
- `DivStats.score()` is non-negative for any tree
- `DivStats + DivStats` is commutative and associative

## Related issues

- Proposal: [kb/proposals/reroot-generic-scoring-architecture.md](../proposals/reroot-generic-scoring-architecture.md)
- Depends on: [kb/tickets/reroot-create-module-with-root-stats-trait.md](reroot-create-module-with-root-stats-trait.md)
- Enables: [kb/tickets/optimize-wire-reroot-with-two-phase-bl-optimization.md](optimize-wire-reroot-with-two-phase-bl-optimization.md)
