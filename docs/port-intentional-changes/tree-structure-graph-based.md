# Graph-based tree structure

This document describes an intentional architectural change affecting tree storage, traversal, node/edge data access, and topology manipulation across all commands.

v1 replaces BioPython's tree representation with a generic directed graph defined in `Graph<N, E, D>` (`#Graph`) in `packages/treetime-graph/src/graph.rs:32-43:`. v0 wraps BioPython's `Phylo.BaseTree.Tree` inside `TreeAnc` (`#TreeAnc`) in `packages/legacy/treetime/treetime/treeanc.py:50-514:`.

## Background: phylogenetic tree representations

A phylogenetic tree is a graph-theoretical representation of evolutionary relationships among taxa [1]. Nodes represent either observed samples (leaves) or inferred ancestors (internal nodes). Edges represent evolutionary lineages with associated branch lengths measured in expected substitutions per site.

Most phylogenetic software uses the Newick format for tree serialization: `(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);`. The Newick Standard was adopted June 26, 1986 by an informal committee at the Society for the Study of Evolution meetings, named after Newick's restaurant in Dover, New Hampshire where they held their final meeting [2]. The format traces mathematically to Arthur Cayley's 1857 observation of the correspondence between trees and nested parentheses [2].

The Newick format encodes parent-child relationships implicitly through nesting. The standard data structure mirrors this: each node stores a list of children, and branch length is stored on the child node [2]. BioPython's `Bio.Phylo.BaseTree.Clade` follows this convention [3].

This child-linked representation has two limitations for algorithms that traverse toward the root. First, accessing a node's parent requires either a separate lookup or an explicit back-pointer. Second, branch data conceptually belongs to the connection between nodes, not to either node individually, but the format forces it onto the child.

## v0: BioPython Clade with monkey-patched attributes

v0 represents trees using BioPython's `Bio.Phylo` module, a unified toolkit for processing, analyzing, and visualizing phylogenetic trees [3]. The tree is a `Phylo.BaseTree.Tree` object. Each node is a `Bio.Phylo.BaseTree.Clade` instance storing children in a `.clades` list attribute [4].

BioPython's data model provides no explicit parent attribute [4]. To find a parent, one must traverse from the root using `get_path()`. v0 adds parent pointers by monkey-patching: `_prepare_nodes()` (`#_prepare_nodes`) in `packages/legacy/treetime/treetime/treeanc.py:461-493:` walks the tree in preorder and sets `c.up = clade` on each child. Every node also receives a back-reference to the owning `TreeAnc` instance via `c.tt = self`.

Algorithm-specific data is attached directly to Clade objects as dynamic attributes. Three properties are patched onto the `Clade` class at module level (`packages/legacy/treetime/treetime/treeanc.py:45-47:`):

- `Clade.sequence` - full sequence array
- `Clade.cseq` - compressed sequence
- `Clade.mutations` - mutations relative to parent

Other attributes are set during algorithm passes without declaration: `node.branch_length_interpolator`, `node.date_constraint`, `node.time_before_present`, `node.marginal_subtree_LH`, `node.joint_Lx`, `node.bad_branch`. Missing attributes cause runtime `AttributeError`.

Traversal uses BioPython methods [4]:

- `tree.find_clades(order='postorder')` for leaves-to-root (depth-first, children before parent)
- `tree.find_clades(order='preorder')` for root-to-leaves (depth-first, parent before children)
- `tree.get_terminals()` for leaf nodes
- `tree.get_nonterminals()` for internal nodes

All traversals are single-threaded and sequential.

## v1: generic directed graph with typed payloads

v1 uses `Graph<N, E, D>` defined in the `treetime-graph` crate. The struct is parameterized over three types:

- `N: GraphNode` - node payload containing algorithm-specific data
- `E: GraphEdge` - edge payload containing branch data
- `D: Sync + Send` - graph-level shared data

Nodes are stored in `Vec<Option<Arc<RwLock<Node<N>>>>>` indexed by `GraphNodeKey(usize)`. Each `Node<N>` (`#Node`) in `packages/treetime-graph/src/node.rs:73-79:` tracks both inbound and outbound edge keys:

```rust
pub struct Node<N: GraphNode> {
  key: GraphNodeKey,
  data: Arc<RwLock<N>>,
  outbound_edges: Vec<GraphEdgeKey>,
  inbound_edges: Vec<GraphEdgeKey>,
  is_visited: AtomicBool,
}
```

Parent access uses `inbound_edges`: `graph.parents_of(node)` (`#parents_of`) in `packages/treetime-graph/src/graph.rs:86-97:` retrieves parent nodes by following inbound edge keys to their source nodes. Child access uses `outbound_edges` similarly via `graph.children_of(node)` (`#children_of`) in `packages/treetime-graph/src/graph.rs:143-154:`.

Edges are stored in `Vec<Option<Arc<RwLock<Edge<E>>>>>`. Each `Edge<E>` (`#Edge`) in `packages/treetime-graph/src/edge.rs:71-83:` stores source and target node keys plus a payload:

```rust
pub struct Edge<E: GraphEdge> {
  key: GraphEdgeKey,
  source: GraphNodeKey,
  target: GraphNodeKey,
  data: Arc<RwLock<E>>,
}
```

Node and edge payloads are concrete types per algorithm stage. `NodeAncestral` (`#NodeAncestral`) and `EdgeAncestral` (`#EdgeAncestral`) in `packages/treetime/src/representation/payload/ancestral.rs:14-72:` provide minimal fields for ancestral reconstruction. `NodeTimetree` (`#NodeTimetree`) and `EdgeTimetree` (`#EdgeTimetree`) in `packages/treetime/src/representation/payload/timetree.rs:18-158:` extend these with time-related fields.

Traits define typed access to domain-specific fields:

- `Named`, `Described` - node identification (`packages/treetime-graph/src/node.rs:11-20:`)
- `Divergence` - evolutionary distance from root (`packages/treetime-graph/src/node.rs:32-35:`)
- `TimeConstraint<T>` - date constraints and bad branch flag (`packages/treetime-graph/src/node.rs:44-49:`)
- `HasBranchLength` - branch length on edges (`packages/treetime-graph/src/edge.rs:13-16:`)
- `ClockMessages<T>` - message-passing fields for clock inference (`packages/treetime-graph/src/edge.rs:20-27:`)

Functions requiring specific data declare trait bounds. Missing data is a compile error.

Traversal methods are built into `Graph`:

- `par_iter_breadth_first_forward` / `par_iter_breadth_first_backward` (`packages/treetime-graph/src/graph_traverse.rs:221-250:`) - parallel breadth-first using rayon, processing independent nodes concurrently within each frontier
- `iter_depth_first_preorder_forward` / `iter_depth_first_postorder_forward` (`packages/treetime-graph/src/graph_traverse.rs:256-321:`) - sequential depth-first
- `iter_breadth_first_forward` / `iter_breadth_first_reverse` (`packages/treetime-graph/src/graph_traverse.rs:327-354:`) - sequential breadth-first

Traversal callbacks receive `GraphNodeForward` (`#GraphNodeForward`) or `GraphNodeBackward` (`#GraphNodeBackward`) in `packages/treetime-graph/src/graph_traverse.rs:18-167:` providing mutable payload access together with read access to parent/child payloads.

## Rationale

BioPython's Clade model has limitations that conflict with v1's requirements:

**No typed payloads.** Algorithm data is stored as dynamic attributes on Clade objects. Any attribute can be set on any node at any time. Missing attributes cause runtime `AttributeError`. v1 uses Rust's type system to ensure each algorithm stage has exactly the data it needs at compile time.

**No explicit edges.** In Bio.Phylo, branch data (length, mutations) is stored on the child node [4]. There is no edge entity. v1's `Edge<E>` provides a natural place for branch-specific data like length distributions and message-passing fields.

**No parent pointers.** BioPython trees are child-linked only [4]. v0 patches parent references manually in `_prepare_nodes()`. v1 stores bidirectional edge references (inbound/outbound) as part of the node structure, making parent access a built-in O(1) operation via `node.inbound()` (`#inbound`) in `packages/treetime-graph/src/node.rs:190-192:`.

**Single-threaded traversal.** BioPython's `find_clades` yields nodes sequentially. v1's parallel breadth-first traversal processes independent nodes concurrently, using frontier-based synchronization to respect data dependencies between tree levels.

**Tree-only topology.** Bio.Phylo assumes strictly tree-shaped data. v1's directed graph supports multiple roots and multiple parents per node, enabling representation of phylogenetic networks. Phylogenetic networks generalize trees to represent reticulate evolution (recombination, hybridization, lateral gene transfer) [5]. The `one_parent_of()` (`#one_parent_of`) method in `packages/treetime-graph/src/graph.rs:121-137:` enforces tree topology when algorithms require it, returning an error for network nodes.

## Practical impact

- Node data access is type-checked at compile time. Algorithms that require specific fields (e.g., time distributions) declare trait bounds, and missing data is a compile error rather than a runtime `AttributeError`.
- Branch data lives on edges. Functions that operate on branches receive `Edge<E>` payloads directly instead of reading attributes from child nodes.
- Traversal order is explicit in the method name. There is no string parameter (`order='postorder'`) that can be misspelled.
- Parallel traversal reduces wall-clock time for large trees. The BFS frontier approach processes all nodes at the same tree depth concurrently.
- The graph structure accepts phylogenetic networks. Trees are a special case where each node has at most one parent.

## References

[1] Qi, J. & Schicho, J. (2020). Five equivalent representations of a phylogenetic tree. arXiv:2011.11774. Provides axiomatic foundation for equivalent tree representations.

[2] Felsenstein, J. The Newick tree format. PHYLIP documentation. https://phylipweb.github.io/phylip/newicktree.html. Documents the Newick Standard adopted June 26, 1986.

[3] Talevich, E., Invergo, B.M., Cock, P.J. & Chapman, B.A. (2012). Bio.Phylo: A unified toolkit for processing, analyzing and visualizing phylogenetic trees in Biopython. BMC Bioinformatics 13:209.

[4] BioPython Bio.Phylo.BaseTree API documentation. https://biopython.org/docs/latest/api/Bio.Phylo.BaseTree.html. Documents the Clade class and traversal methods.

[5] Cardona, G., Rossello, F. & Valiente, G. (2008). Extended Newick: it is time for a standard representation of phylogenetic networks. BMC Bioinformatics 9:532. PMC2621367.
