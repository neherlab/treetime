# Graph-based tree structure

| Property    | Value                                                                                                                                                                                                                                                          |
| ----------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Intentional architectural change                                                                                                                                                                                                                               |
| v1 Location | `Graph<N, E, D>` (`#Graph`) in [`packages/treetime-graph/src/graph.rs#L32-L43`](../../packages/treetime-graph/src/graph.rs#L32-L43), `Node<N>` (`#Node`) in [`packages/treetime-graph/src/node.rs#L73-L79`](../../packages/treetime-graph/src/node.rs#L73-L79) |
| v0 Location | `TreeAnc` (`#TreeAnc`) in [`packages/legacy/treetime/treetime/treeanc.py#L50-L514`](../../packages/legacy/treetime/treetime/treeanc.py#L50-L514), wrapping `Bio.Phylo.BaseTree.Tree`                                                                           |
| Affects     | Tree storage, traversal, node/edge data access, topology manipulation                                                                                                                                                                                          |
| Commands    | All commands                                                                                                                                                                                                                                                   |

## What v0 does

v0 represents phylogenetic trees using Biopython's `Bio.Phylo` module. The tree is a `Phylo.BaseTree.Tree` object, and each node is a `Bio.Phylo.BaseTree.Clade` instance.

Clade objects store children in a `.clades` list attribute. Parent access does not exist in Biopython's data model, so v0 adds it by monkey-patching: `_prepare_nodes()` walks the tree in preorder and sets `c.up = clade` on each child. Every node also receives a back-reference to the owning `TreeAnc` instance via `c.tt = self`.

Algorithm-specific data (compressed sequences, mutation lists, date constraints, time distributions) is attached directly to Clade objects as dynamic attributes. For example, `Clade.cseq` and `Clade.mutations` are added as Python properties patched onto the `Clade` class at module level. Other attributes like `node.branch_length_interpolator`, `node.date_constraint`, `node.time_before_present`, and `node.bad_branch` are set during various algorithm passes.

Traversal uses Biopython methods:

- `tree.find_clades(order='postorder')` for leaves-to-root
- `tree.find_clades(order='preorder')` for root-to-leaves
- `tree.get_terminals()` for leaf nodes
- `tree.get_nonterminals()` for internal nodes

All traversals are single-threaded and sequential.

## What v1 does

v1 uses a generic directed graph `Graph<N, E, D>` defined in the `treetime-graph` crate. The struct is parameterized over three types:

- `N: GraphNode` - node payload (the algorithm-specific data attached to each node)
- `E: GraphEdge` - edge payload (branch lengths, message-passing distributions)
- `D: Sync + Send` - graph-level shared data

Nodes and edges are stored in `Vec<Option<Arc<RwLock<Node<N>>>>>` and `Vec<Option<Arc<RwLock<Edge<E>>>>>`, indexed by `GraphNodeKey(usize)` and `GraphEdgeKey(usize)`. Each `Node<N>` tracks both inbound and outbound edge keys, providing O(1) parent and child access without monkey-patching.

Node and edge payloads are separate types defined per algorithm stage. Traits like `Named`, `Divergence`, `TimeConstraint`, `HasBranchLength`, and `ClockMessages` define typed access to domain-specific fields. Different commands use different payload types: `NodeAncestral`/`EdgeAncestral` for ancestral reconstruction, `NodeTimetree`/`EdgeTimetree` for timetree inference.

Traversal methods are built into the graph:

- `par_iter_breadth_first_forward` / `par_iter_breadth_first_backward` - parallel breadth-first using rayon, processing independent nodes concurrently within each frontier
- `iter_depth_first_preorder_forward` / `iter_depth_first_postorder_forward` - sequential depth-first
- `iter_breadth_first_forward` / `iter_breadth_first_reverse` - sequential breadth-first

Traversal callbacks receive `GraphNodeForward` or `GraphNodeBackward` wrapper structs that provide mutable payload access together with read access to parent/child payloads, enforced at compile time.

## Why v1 changes this

Biopython's Clade model has limitations that conflict with v1's design:

- **No typed payloads.** Algorithm data is stored as dynamic attributes on Clade objects. Any attribute can be set on any node at any time, and missing attributes cause runtime `AttributeError`. v1 uses Rust's type system to ensure each algorithm stage has exactly the data it needs.

- **No explicit edges.** In Bio.Phylo, branch data (length, mutations) is stored on the child node. There is no edge entity. v1's `Edge<E>` struct provides a natural place for branch-specific data like length distributions and message-passing fields.

- **No parent pointers.** Biopython trees are child-linked only. v0 patches parent references manually. v1 stores bidirectional edge references (inbound/outbound) as part of the node structure, making parent access a built-in operation.

- **Single-threaded traversal.** Biopython's `find_clades` yields nodes sequentially. v1's parallel breadth-first traversal processes independent nodes concurrently, using frontier-based synchronization to respect data dependencies.

- **Tree-only topology.** Bio.Phylo assumes strictly tree-shaped data. v1's directed graph supports multiple roots and multiple parents per node, allowing representation of phylogenetic networks (ARGs) with the same data structure.

## Practical impact

- Node data access is type-checked at compile time. Algorithms that require specific fields (e.g. time distributions) declare trait bounds, and missing data is a compile error rather than a runtime `AttributeError`.
- Branch data lives on edges. Functions that operate on branches receive `Edge<E>` payloads directly instead of reading attributes from child nodes.
- Traversal order is explicit in the method name. There is no string parameter (`order='postorder'`) that can be misspelled or confused.
- Parallel traversal reduces wall-clock time for large trees. The BFS frontier approach processes all nodes at the same tree depth concurrently.
- The graph structure accepts phylogenetic networks. Trees are a special case where each node has at most one parent. The `one_parent_of()` method enforces tree topology when algorithms require it, returning an error for network nodes.
