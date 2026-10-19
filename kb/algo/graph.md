# Graph Traversal Algorithms

## Parallel Tree Scans

Two parallel scans on a `TreeView` run every pass that follows the tree topology (<a id="cite-1"></a>[Leiserson and Schardl 2010](https://doi.org/10.1145/1810479.1810534) [[1](#ref-1)]): the Fitch, dense, discrete, and sparse marginal, clock regression, and timetree passes. `scan_down` visits each node after its parent (forward passes). `scan_up` visits each node after all of its children (backward passes).

Each scan takes one input per node by value. The step of a node returns the output of the node and one message: `scan_down` sends the message to every child, and `scan_up` sends it to the parent. When a pass needs the same value as output and as message, it shares the value through `Arc` and unwraps it after the scan, so the scan does not clone sequence or probability-matrix state. `join_parent_edges()` pairs each node input with the input of its parent edge, and `split_parent_edges()` splits the outputs into a node map and a parent-edge map. Each result map is bulk-built in one step from its entries in key order, which fills the map nodes; inserts in key order would leave them about half full, and the allocator would keep more freed memory.

The task that makes a node ready runs that node, so no task waits on a queue or on a per-level barrier:

- `scan_down` runs in one `rayon::scope`. A task visits its node, spawns one task per child except the first, and continues with the first child
- `scan_up` starts one task per worker thread, and no more tasks than leaves. Each task takes the next leaf in preorder from a shared counter over the leaf list and walks up from it: it puts its message into the slot of its child position in the inbox of the parent, the delivery that fills the last slot continues with the parent, and every other delivery ends the walk, so the task takes the next leaf
- `scan_up` gives the step the messages of the children in the order of `TreeView::children()`, never in arrival order, so floating-point products and sums are the same for every thread count
- The first error wins. Later visits stop, and in `scan_up` the ancestors of a failed node do not run
- Both scans loop instead of recursing, and no task waits for another task, so stack use depends neither on the tree height nor on work stealing

v1: [`packages/treetime-graph/src/tree_scan.rs`](../../packages/treetime-graph/src/tree_scan.rs), tree shape in [`packages/treetime-graph/src/tree_view.rs`](../../packages/treetime-graph/src/tree_view.rs).

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
| [`packages/treetime-graph/src/tree_scan.rs`](../../packages/treetime-graph/src/tree_scan.rs)               | Parallel down and up tree scans    |
| [`packages/treetime-graph/src/tree_view.rs`](../../packages/treetime-graph/src/tree_view.rs)               | Tree shape checks and child order  |
| [`packages/treetime-graph/src/graph_traverse.rs`](../../packages/treetime-graph/src/graph_traverse.rs)     | DFS preorder/postorder, serial BFS |
| [`packages/treetime-graph/src/reachability.rs`](../../packages/treetime-graph/src/reachability.rs)         | Forward reachability query         |
| [`packages/treetime-graph/src/graph_ops.rs`](../../packages/treetime-graph/src/graph_ops.rs)               | Edge collapse                      |
| [`packages/treetime-graph/src/topology_order.rs`](../../packages/treetime-graph/src/topology_order.rs)     | Topology ordering                  |
