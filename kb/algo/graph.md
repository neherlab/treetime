# Graph Traversal Algorithms

## Graph Schedule

`GraphSchedule::new()` indexes a graph once and gives every pass the same validated order. It works on rooted directed acyclic graphs (DAGs): a node may have several parents, and the graph may have several roots. The schedule borrows the graph, so the topology cannot change while a schedule exists.

- **Adjacency**: the parents of each node in `inbound()` order and the children in `outbound()` order, each with its edge key, in flat arrays indexed by node key. `parent_slot(edge)` is the position of an edge among the parents of its child, `child_slot(edge)` its position among the children of its parent. The slots belong to edges, because a shared child has a different position in each parent
- **Validation**: the graph must have a root ("The graph has no root: every node has a parent"), and the preorder must visit every node, which holds exactly when the graph has no cycle. Otherwise the error names a node on a cycle ("The graph has a cycle through node {key}"), found by a walk from an unvisited node along unvisited parents until a node repeats
- **Preorder**: a stack starts with the roots in `graph.roots` order; a node is pushed when its last parent is visited, so every node comes after all its parents. On a tree this is the depth-first preorder with children in `outbound()` order
- **Postorder**: the depth-first finish order (<a id="cite-1"></a>[Cormen et al. 2022](https://mitpress.mit.edu/9780262046305/introduction-to-algorithms/) [[1](#ref-1)], Chapter 22.3) with explicit enter and exit frames on one stack. Each node appears once, after all its children. Reversed preorder is a different order on trees (siblings come in reverse), so the schedule keeps both
- **Sequential loops**: every serial pass iterates `schedule.preorder()` or `schedule.postorder()` and reads neighbors through `parents(key)` and `children(key)`, which are slices into the schedule, so a visit allocates nothing
- **Tree check**: `require_tree()` checks one root and at most one parent per node. Tree-only algorithms (Fitch, marginal, clock, timetree) call it once at their start
- **Tree handle**: `into_tree()` runs the check and returns a `TreeSchedule` with `root()`, `parent()`, `children()`, and `preorder()`. Only the output formats that can hold only trees (Newick, Nexus, Auspice JSON, UShER MAT, augur node data, per-node tables) take it, through `AnnotatedTreeView`. The `forbidden_item` lint of `dev/lints/dylint-custom` (configured in `dylint.toml`) rejects `TreeSchedule` and `into_tree()` in production code of the `treetime` crate, so core algorithms stay on the graph schedule

v1: [`packages/treetime-graph/src/schedule.rs`](../../packages/treetime-graph/src/schedule.rs).

---

## Parallel Scans

Two parallel scans on a `GraphSchedule` run every pass that follows the topology (<a id="cite-2"></a>[Leiserson and Schardl 2010](https://doi.org/10.1145/1810479.1810534) [[2](#ref-2)]): the Fitch, dense, discrete, and sparse marginal, clock regression, and timetree passes. `scan_down` visits each node after all its parents (forward passes). `scan_up` visits each node after all its children (backward passes).

Each scan takes one input per node by value. The step of a node gets its input and a `Vec` of messages, one per parent in `parents(key)` order (`scan_down`, empty at a root) or one per child in `children(key)` order (`scan_up`, empty at a leaf), and returns the output of the node and one message. The message goes to every child (`scan_down`) or every parent (`scan_up`), a clone per edge. A tree pass takes its single parent message with `messages.pop()`. When a pass needs the same value as output and as message, it shares the value through `Arc` and unwraps it after the scan (`unwrap_shared()`), so the scan does not clone sequence or probability-matrix state. `join_parent_edges()` pairs each node input with one edge input per parent, and `split_parent_edges()` splits the outputs into a node map and a parent-edge map and checks one edge value per parent edge. Each result map is bulk-built in one step from its entries in key order, which fills the map nodes; inserts in key order would leave them about half full, and the allocator would keep more freed memory.

The task that makes a node ready runs that node, so no task waits on a queue or on a per-level barrier:

- Each node has an inbox with one slot per input edge (`parent_slot` for `scan_down`, `child_slot` for `scan_up`). The delivery that fills the last slot continues with that node; a node with one input skips the inbox
- `scan_down` runs in one `rayon::scope` and starts one task per root. A task visits its node, continues with the first child that becomes ready, and spawns one task per further ready child
- `scan_up` starts one task per worker thread, and no more tasks than leaves. Each task takes the next leaf in preorder from a shared counter over the leaf list and walks up from it while its deliveries complete an inbox; then it takes the next leaf
- The step gets its messages in edge order, never in arrival order, so floating-point products and sums are the same for every thread count
- The first error wins. Later visits stop, and nodes that wait for a failed node do not run
- Both scans loop instead of recursing, and no task waits for another task, so stack use depends neither on the graph depth nor on work stealing. A step runs no rayon parallelism of its own, because a waiting `join` inside a step runs stolen work on top of its stack

v1: [`packages/treetime-graph/src/schedule_scan.rs`](../../packages/treetime-graph/src/schedule_scan.rs).

---

## Reachability

Serial forward-reachability query: whether a directed path runs from one node to another along edge directions. A depth-first walk over outbound edges with a visited set for cycle termination.

v1: [`packages/treetime-graph/src/reachability.rs`](../../packages/treetime-graph/src/reachability.rs).

---

## Path Finding

O(h) parent-chain traversal for finding paths between nodes in the tree, where h is the tree height. Used by rerooting to trace the path from the current root to the new root position.

v1: [`packages/treetime-graph/src/graph.rs#L161-L211`](../../packages/treetime-graph/src/graph.rs#L161-L211).

---

## Edge Contraction

Merges a target node into its source (parent) node (<a id="cite-3"></a>[Diestel 2017](https://doi.org/10.1007/978-3-662-53622-3) [[3](#ref-3)], Chapter 1) by reparenting all of the target's children to the source and removing the contracted edge. Used to clean up single-child internal nodes after polytomy resolution and rerooting.

v1: [`packages/treetime-graph/src/graph_ops.rs#L139-L186`](../../packages/treetime-graph/src/graph_ops.rs#L139-L186).

---

## Topology Ordering

Output-boundary topology ordering computes a deterministic child order before writing graph-backed tree files. The default order matches TreeTime v0 ladderization by sorting siblings by reachable leaf count in ascending order. Additional presets preserve input order, reverse ladderization, sort by subtree height, sort by leaf labels, or sort against an explicit target leaf order. Target-order scores use mean or median target position across each subtree's descendant leaves. DAG inputs count each reachable leaf once per child and cyclic inputs fail before output with the errors of the graph schedule, whose postorder gives the order of the key computation.

v1: [`packages/treetime-graph/src/topology_order.rs`](../../packages/treetime-graph/src/topology_order.rs).

---

## References

- <a id="ref-1"></a>Cormen, Thomas H., Charles E. Leiserson, Ronald L. Rivest, and Clifford Stein. 2022. _Introduction to Algorithms._ 4th ed. MIT Press. ISBN 978-0-262-04630-5. https://mitpress.mit.edu/9780262046305/introduction-to-algorithms/ [↩](#cite-1)
- <a id="ref-2"></a>Leiserson, Charles E., and Tao B. Schardl. 2010. "A Work-Efficient Parallel Breadth-First Search Algorithm (or How to Cope with the Nondeterminism of Reducers)." In _Proceedings of the 22nd ACM Symposium on Parallelism in Algorithms and Architectures (SPAA),_ 303-314. https://doi.org/10.1145/1810479.1810534 [↩](#cite-2)
- <a id="ref-3"></a>Diestel, Reinhard. 2017. _Graph Theory._ 5th ed. Springer. ISBN 978-3-662-53621-6. https://doi.org/10.1007/978-3-662-53622-3 [↩](#cite-3)

---

## File Index

| File                                                                                                   | Algorithms                                  |
| ------------------------------------------------------------------------------------------------------ | ------------------------------------------- |
| [`packages/treetime-graph/src/schedule.rs`](../../packages/treetime-graph/src/schedule.rs)             | Graph schedule: orders, checks, tree handle |
| [`packages/treetime-graph/src/schedule_scan.rs`](../../packages/treetime-graph/src/schedule_scan.rs)   | Parallel down and up scans                  |
| [`packages/treetime-graph/src/reachability.rs`](../../packages/treetime-graph/src/reachability.rs)     | Forward reachability query                  |
| [`packages/treetime-graph/src/graph_ops.rs`](../../packages/treetime-graph/src/graph_ops.rs)           | Edge collapse                               |
| [`packages/treetime-graph/src/topology_order.rs`](../../packages/treetime-graph/src/topology_order.rs) | Topology ordering                           |
