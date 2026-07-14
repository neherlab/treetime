# Dependency-frontier schedulers are duplicated

Graph traversal and indexed partition passes each implement dependency counts, ready-frontier advancement, parallel execution, and cancellation behavior. Their ordering and failure semantics can diverge as either implementation changes.

The graph layer builds and executes frontiers in [`packages/treetime-graph/src/graph_traverse.rs`](../../packages/treetime-graph/src/graph_traverse.rs), while indexed inference owns a second frontier representation and traversal loop in [`packages/treetime/src/partition/indexed_pass.rs`](../../packages/treetime/src/partition/indexed_pass.rs).

## Potential solutions

- O1. Extract a generic dependency scheduler with domain callbacks.
- O2. Choose one existing scheduler and adapt the other domain directly to it. This is preferable if its current API already expresses all failure/commit invariants.

## Recommendation

Extract one generic dependency scheduler. Domain layers provide the dependency graph, work callback, and deterministic commit callback; the scheduler owns readiness, cancellation, and completion accounting.

## Related issues

- [N-graph-traverse-quadratic-iter-children-arc.md](N-graph-traverse-quadratic-iter-children-arc.md)
