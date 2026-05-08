# Graph Traversal Algorithms

[Back to index](README.md)

## Parallel BFS

Frontier-based parallel breadth-first traversal using Rayon for work distribution across CPU cores. Each BFS level is processed in parallel, with synchronization between levels. Used for Fitch parsimony and marginal reconstruction where all nodes at the same depth can be processed independently.

v1: [`packages/treetime-graph/src/breadth_first.rs#L139-L202`](../../packages/treetime-graph/src/breadth_first.rs#L139-L202).

Reference: Leiserson & Schardl (2010). "A Work-Efficient Parallel BFS Algorithm." SPAA 2010. doi:10.1145/1810479.1810534

---

## DFS Preorder/Postorder

Iterative stack-based depth-first traversal supporting both preorder (parent before children) and postorder (children before parent) visitation. Postorder is used for leaf-to-root passes (Fitch backward, marginal backward, clock regression backward). Preorder is used for root-to-leaf passes (Fitch forward, marginal forward, clock regression forward).

v1: [`packages/treetime-graph/src/graph_traverse.rs#L256-L321`](../../packages/treetime-graph/src/graph_traverse.rs#L256-L321).

Reference: Cormen, Leiserson, Rivest & Stein (2022). "Introduction to Algorithms." 4th ed. MIT Press, Chapter 22.3. ISBN 978-0-262-04630-5

---

## Path Finding

O(h) parent-chain traversal for finding paths between nodes in the tree, where h is the tree height. Used by rerooting to trace the path from the current root to the new root position.

v1: [`packages/treetime-graph/src/graph.rs#L316-L370`](../../packages/treetime-graph/src/graph.rs#L316-L370).

---

## Edge Contraction

Merges a target node into its source (parent) node by reparenting all of the target's children to the source and removing the contracted edge. Used to clean up single-child internal nodes after polytomy resolution and rerooting.

v1: [`packages/treetime-graph/src/graph_ops.rs#L146-L206`](../../packages/treetime-graph/src/graph_ops.rs#L146-L206).

Reference: Diestel (2017). "Graph Theory." 5th ed. Springer, Chapter 1. ISBN 978-3-662-53621-6

---

## File Index

| File                                                                                                   | Algorithms                             |
| ------------------------------------------------------------------------------------------------------ | -------------------------------------- |
| [`packages/treetime-graph/src/breadth_first.rs`](../../packages/treetime-graph/src/breadth_first.rs)   | Parallel BFS                           |
| [`packages/treetime-graph/src/graph_traverse.rs`](../../packages/treetime-graph/src/graph_traverse.rs) | DFS preorder/postorder, sequential BFS |
| [`packages/treetime-graph/src/find_paths.rs`](../../packages/treetime-graph/src/find_paths.rs)         | Path finding                           |
| [`packages/treetime-graph/src/graph_ops.rs`](../../packages/treetime-graph/src/graph_ops.rs)           | Edge collapse                          |
