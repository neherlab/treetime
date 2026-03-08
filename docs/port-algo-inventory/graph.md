# Graph Traversal Algorithms

[Back to index](../index.md)

## Parallel BFS

| Property    | Value                                                                                                                       |
| ----------- | --------------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known                                                                                                                  |
| v1 Location | [`packages/treetime-graph/src/breadth_first.rs#L139-L202`](../../../packages/treetime-graph/src/breadth_first.rs#L139-L202) |
| Reference   | Leiserson & Schardl (2010). "A Work-Efficient Parallel BFS Algorithm." SPAA 2010                                            |

Frontier-based parallel traversal using Rayon.

---

## DFS Preorder/Postorder

| Property    | Value                                                                                                                         |
| ----------- | ----------------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known                                                                                                                    |
| v1 Location | [`packages/treetime-graph/src/graph_traverse.rs#L256-L321`](../../../packages/treetime-graph/src/graph_traverse.rs#L256-L321) |
| Reference   | CLRS Chapter 22.3                                                                                                             |

Iterative stack-based traversal.

---

## Path Finding

| Property    | Value                                                                                                       |
| ----------- | ----------------------------------------------------------------------------------------------------------- |
| Type        | Standard tree operation                                                                                     |
| v1 Location | [`packages/treetime-graph/src/graph.rs#L316-L370`](../../../packages/treetime-graph/src/graph.rs#L316-L370) |

O(h) parent-chain traversal.

---

## Edge Contraction

| Property    | Value                                                                                                               |
| ----------- | ------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known                                                                                                          |
| v1 Location | [`packages/treetime-graph/src/graph_ops.rs#L146-L206`](../../../packages/treetime-graph/src/graph_ops.rs#L146-L206) |
| Reference   | Diestel, "Graph Theory." Chapter 1                                                                                  |

Merges target node into source node.

---

## File Index

| File                                                                                                      | Algorithms                             |
| --------------------------------------------------------------------------------------------------------- | -------------------------------------- |
| [`packages/treetime-graph/src/breadth_first.rs`](../../../packages/treetime-graph/src/breadth_first.rs)   | Parallel BFS                           |
| [`packages/treetime-graph/src/graph_traverse.rs`](../../../packages/treetime-graph/src/graph_traverse.rs) | DFS preorder/postorder, sequential BFS |
| [`packages/treetime-graph/src/find_paths.rs`](../../../packages/treetime-graph/src/find_paths.rs)         | Path finding                           |
| [`packages/treetime-graph/src/graph_ops.rs`](../../../packages/treetime-graph/src/graph_ops.rs)           | Edge collapse                          |
