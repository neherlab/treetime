# treetime-graph

Generic directed graph data structure for phylogenetic trees. Provides thread-safe nodes and edges with user-defined payloads, plus synchronous and parallel traversal algorithms.

## Key types

| Type             | Description                                                                                  |
| ---------------- | -------------------------------------------------------------------------------------------- |
| `Graph<N, E, D>` | Directed graph parameterized by node payload `N`, edge payload `E`, and graph-level data `D` |
| `Node<N>`        | Graph node wrapping a payload `N` with inbound/outbound edge tracking                        |
| `Edge<E>`        | Directed edge connecting source to target node, wrapping a payload `E`                       |
| `GraphNodeKey`   | Newtype index into the node storage (`usize`)                                                |
| `GraphEdgeKey`   | Newtype index into the edge storage (`usize`)                                                |

All nodes and edges are stored as `Arc<RwLock<_>>` (using `parking_lot`) for concurrent access. Type aliases `SafeNode`, `SafeEdge`, `SafeNodeRef`, etc. wrap the lock guard types.

## Node traits

- `GraphNode` - marker trait (`Clone + Debug + Sync + Send`)
- `Named` - read/write node name
- `Described` - read/write node description
- `Divergence` - evolutionary distance from root
- `Outlier` - outlier flag
- `TimeConstraint<T>` - time distribution and bad-branch flag

Composite traits (`NodeAncestralOps`, `NodeOptimizeOps`) combine these for specific algorithms.

## Edge traits

- `GraphEdge` - marker trait (`Clone + Debug + Sync + Send`)
- `HasBranchLength` - read/write branch length
- `ClockMessages<T>` - message passing fields for clock inference
- `BranchDistribution<T>` - branch length distribution and parent message
- `TimeLength` - time-scaled branch length

## Graph operations

Build a graph by adding nodes and edges, then call `build()` to compute root and leaf sets:

```rust
let mut graph = Graph::<MyNode, MyEdge>::new();
let a = graph.add_node(MyNode::new("A"));
let b = graph.add_node(MyNode::new("B"));
graph.add_edge(a, b, MyEdge::default())?;
graph.build()?;
```

Mutation operations: `add_node`, `add_edge`, `remove_node`, `remove_edge`, `collapse_edge`.

Query operations: `get_node`, `get_edge`, `find_node`, `parents_of`, `children_of`, `get_roots`, `get_leaves`, `get_internal_nodes`, `path_from_root_to_node`.

## Traversal

### Synchronous (single-threaded)

| Method                               | Direction      | Order          |
| ------------------------------------ | -------------- | -------------- |
| `iter_depth_first_preorder_forward`  | root to leaves | DFS pre-order  |
| `iter_depth_first_postorder_forward` | leaves to root | DFS post-order |
| `iter_breadth_first_forward`         | root to leaves | BFS            |
| `iter_breadth_first_reverse`         | leaves to root | BFS reverse    |

Each traversal method takes a closure receiving `GraphNodeForward` or `GraphNodeBackward`, which provides mutable access to the current node's payload plus read access to parent/child payloads and edges.

### Parallel (rayon-based)

| Method                            | Direction      |
| --------------------------------- | -------------- |
| `par_iter_breadth_first_forward`  | root to leaves |
| `par_iter_breadth_first_backward` | leaves to root |

Parallel traversal processes each frontier (set of nodes whose dependencies are resolved) concurrently using rayon. Returns `GraphTraversalContinuation` to allow early termination.

## Path finding

- `find_paths` - find all edges on paths between two nodes
- `exists_forward_path_between` / `exists_backward_path_between` - check path existence
- `path_from_root_to_node` / `path_from_node_to_node` - collect nodes along a path

## Utilities

- `assign_node_names` - assign auto-generated names (`NODE_0000000`, ...) to unnamed nodes
- `invert_edge` - reverse an edge's direction, updating adjacency lists on both endpoints
