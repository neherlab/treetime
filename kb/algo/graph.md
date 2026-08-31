# Graph Traversal Algorithms

[Back to index](README.md)

## Work-first Parallel Pass

The single parallel tree traversal (<a id="cite-1"></a>[Leiserson and Schardl 2010](https://doi.org/10.1145/1810479.1810534) [[1](#ref-1)]). Fitch, marginal, clock, and timetree passes arrange partition payloads in stable indexed slots and run them in one Rayon scope. Roots seed forward passes and leaves seed backward passes. Each completed node decrements an atomic prerequisite counter for its successors; a successor enters the shared work queue as soon as its counter reaches zero, so a ready node runs without waiting on a per-level barrier.

Workers take ownership of one mutable slot and publish it once after processing. Domain callbacks receive read-only access to completed dependency payloads and the writable current slot. Payloads are moved into and out of the pass, so scheduling does not clone sequence or probability-matrix state, while parent-before-child (forward) or child-before-parent (backward) ordering is preserved.

v1: [`packages/treetime-graph/src/pass.rs`](../../packages/treetime-graph/src/pass.rs), scheduler in [`packages/treetime-graph/src/dependency_queue.rs`](../../packages/treetime-graph/src/dependency_queue.rs).

---

## DFS Preorder/Postorder

Iterative stack-based depth-first traversal (<a id="cite-2"></a>[Cormen et al. 2022](https://mitpress.mit.edu/9780262046305/introduction-to-algorithms/) [[2](#ref-2)], Chapter 22.3) supporting both preorder (parent before children) and postorder (children before parent) visitation. These serial walks are used where the callback must capture mutable outer state (a parallel callback cannot): deterministic-order output, IO serialization, and local aggregation.

v1: [`packages/treetime-graph/src/graph_traverse.rs`](../../packages/treetime-graph/src/graph_traverse.rs).

---

## Reachability

Serial forward-reachability query: whether a directed path runs from one node to another along edge directions. A depth-first walk over outbound edges with a visited set for cycle termination.

v1: [`packages/treetime-graph/src/reachability.rs`](../../packages/treetime-graph/src/reachability.rs).

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
- <a id="ref-2"></a>Cormen, Thomas H., Charles E. Leiserson, Ronald L. Rivest, and Clifford Stein. 2022. _Introduction to Algorithms._ 4th ed. MIT Press. ISBN 978-0-262-04630-5. https://mitpress.mit.edu/9780262046305/introduction-to-algorithms/ [↩](#cite-2)
- <a id="ref-3"></a>Diestel, Reinhard. 2017. _Graph Theory._ 5th ed. Springer. ISBN 978-3-662-53621-6. https://doi.org/10.1007/978-3-662-53622-3 [↩](#cite-3)

---

## File Index

| File                                                                                                       | Algorithms                         |
| ---------------------------------------------------------------------------------------------------------- | ---------------------------------- |
| [`packages/treetime-graph/src/pass.rs`](../../packages/treetime-graph/src/pass.rs)                         | Work-first parallel pass           |
| [`packages/treetime-graph/src/dependency_queue.rs`](../../packages/treetime-graph/src/dependency_queue.rs) | Work-first scheduler               |
| [`packages/treetime-graph/src/graph_traverse.rs`](../../packages/treetime-graph/src/graph_traverse.rs)     | DFS preorder/postorder, serial BFS |
| [`packages/treetime-graph/src/reachability.rs`](../../packages/treetime-graph/src/reachability.rs)         | Forward reachability query         |
| [`packages/treetime-graph/src/graph_ops.rs`](../../packages/treetime-graph/src/graph_ops.rs)               | Edge collapse                      |
| [`packages/treetime-graph/src/topology_order.rs`](../../packages/treetime-graph/src/topology_order.rs)     | Topology ordering                  |
