# Graph traversal allocates on every visit

`GraphNodeForward::new()` and `GraphNodeBackward::new()` collect the parent and child keys of each visited node into new `Vec`s ([packages/treetime-graph/src/graph_traverse.rs#L83-L84](../../packages/treetime-graph/src/graph_traverse.rs#L83-L84), [#L112-L113](../../packages/treetime-graph/src/graph_traverse.rs#L112-L113)), two allocations per node per traversal.

Required behavior: no per-visit allocation; callers that need the neighbor keys read them from the graph.
