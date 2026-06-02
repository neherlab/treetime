# Graph Traversal Algorithms

[Back to index](README.md)

## Parallel BFS

Frontier-based parallel breadth-first traversal (<a id="cite-1"></a>[Leiserson and Schardl 2010](https://doi.org/10.1145/1810479.1810534) [[1](#ref-1)]) using Rayon for work distribution across CPU cores. Each BFS level is processed in parallel, with synchronization between levels. Used for Fitch parsimony and marginal reconstruction where all nodes at the same depth can be processed independently.

v1: [`packages/treetime-graph/src/breadth_first.rs#L139-L202`](../../packages/treetime-graph/src/breadth_first.rs#L139-L202).

---

## DFS Preorder/Postorder

Iterative stack-based depth-first traversal (<a id="cite-2"></a>[Cormen et al. 2022](https://doi.org/10.7551/mitpress/13309.001.0001) [[2](#ref-2)], Chapter 22.3) supporting both preorder (parent before children) and postorder (children before parent) visitation. Postorder is used for leaf-to-root passes (Fitch backward, marginal backward, clock regression backward). Preorder is used for root-to-leaf passes (Fitch forward, marginal forward, clock regression forward).

v1: [`packages/treetime-graph/src/graph_traverse.rs#L256-L321`](../../packages/treetime-graph/src/graph_traverse.rs#L256-L321).

---

## Path Finding

O(h) parent-chain traversal for finding paths between nodes in the tree, where h is the tree height. Used by rerooting to trace the path from the current root to the new root position.

v1: [`packages/treetime-graph/src/graph.rs#L316-L370`](../../packages/treetime-graph/src/graph.rs#L316-L370).

---

## Edge Contraction

Merges a target node into its source (parent) node (<a id="cite-3"></a>[Diestel 2017](https://doi.org/10.1007/978-3-662-53622-3) [[3](#ref-3)], Chapter 1) by reparenting all of the target's children to the source and removing the contracted edge. Used to clean up single-child internal nodes after polytomy resolution and rerooting.

v1: [`packages/treetime-graph/src/graph_ops.rs#L146-L206`](../../packages/treetime-graph/src/graph_ops.rs#L146-L206).

---

## Topology Ordering

Output-boundary topology ordering computes a deterministic child order before writing graph-backed tree files. The default order matches TreeTime v0 ladderization by sorting siblings by reachable leaf count in ascending order. Additional presets preserve input order, reverse ladderization, sort by subtree height, sort by leaf labels, or sort against an explicit target leaf order. Target-order scores use mean or median target position across each subtree's descendant leaves. DAG inputs count each reachable leaf once per child and cyclic inputs fail before output.

v1: [`packages/treetime-graph/src/topology_order.rs`](../../packages/treetime-graph/src/topology_order.rs).

---

## References

- <a id="ref-1"></a>Leiserson, Charles E., and Tao B. Schardl. 2010. "A Work-Efficient Parallel Breadth-First Search Algorithm (or How to Cope with the Nondeterminism of Reducers)." In _Proceedings of the 22nd ACM Symposium on Parallelism in Algorithms and Architectures (SPAA),_ 303-314. https://doi.org/10.1145/1810479.1810534 [↩](#cite-1)
- <a id="ref-2"></a>Cormen, Thomas H., Charles E. Leiserson, Ronald L. Rivest, and Clifford Stein. 2022. _Introduction to Algorithms._ 4th ed. MIT Press. ISBN 978-0-262-04630-5. [↩](#cite-2)
- <a id="ref-3"></a>Diestel, Reinhard. 2017. _Graph Theory._ 5th ed. Springer. ISBN 978-3-662-53621-6. https://doi.org/10.1007/978-3-662-53622-3 [↩](#cite-3)

---

## File Index

| File                                                                                                   | Algorithms                             |
| ------------------------------------------------------------------------------------------------------ | -------------------------------------- |
| [`packages/treetime-graph/src/breadth_first.rs`](../../packages/treetime-graph/src/breadth_first.rs)   | Parallel BFS                           |
| [`packages/treetime-graph/src/graph_traverse.rs`](../../packages/treetime-graph/src/graph_traverse.rs) | DFS preorder/postorder, sequential BFS |
| [`packages/treetime-graph/src/find_paths.rs`](../../packages/treetime-graph/src/find_paths.rs)         | Path finding                           |
| [`packages/treetime-graph/src/graph_ops.rs`](../../packages/treetime-graph/src/graph_ops.rs)           | Edge collapse                          |
| [`packages/treetime-graph/src/topology_order.rs`](../../packages/treetime-graph/src/topology_order.rs) | Topology ordering                      |
