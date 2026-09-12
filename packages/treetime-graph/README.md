# treetime-graph

Directed graph data structure for phylogenetic trees. Provides thread-safe node and edge structure with graph-level command data `D`, plus synchronous and parallel traversal algorithms. Per-node and per-edge data lives in external value maps keyed by `GraphNodeKey`/`GraphEdgeKey`.

## Key types

| Type           | Description                                                         |
| -------------- | ------------------------------------------------------------------- |
| `Graph<D>`     | Directed graph carrying graph-level command data `D` (default `()`) |
| `Node`         | Graph node with inbound/outbound edge tracking                      |
| `Edge`         | Directed edge connecting a source to a target node                  |
| `GraphNodeKey` | Newtype index into the node storage (`usize`)                       |
| `GraphEdgeKey` | Newtype index into the edge storage (`usize`)                       |

All nodes and edges are stored as `Arc<RwLock<_>>` (using `parking_lot`) for concurrent traversal. Type aliases `SafeNode`, `SafeEdge`, `SafeNodeRef`, etc. wrap the lock guard types.

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

Each traversal method takes a closure receiving `GraphNodeForward` or `GraphNodeBackward`, which exposes the current node's key and its parent/child keys and edge keys, so the closure can index the external value maps it operates on.

### Parallel (rayon-based)

| Method                            | Direction      |
| --------------------------------- | -------------- |
| `par_iter_breadth_first_forward`  | root to leaves |
| `par_iter_breadth_first_backward` | leaves to root |

Parallel traversal processes each frontier (set of nodes whose dependencies are resolved) concurrently using rayon. Returns `GraphTraversalContinuation` to allow early termination.

Each traversal invocation tracks completed nodes locally by `GraphNodeKey`. Concurrent traversals over one graph remain independent, and early termination or errors require no graph reset.

## Path finding

- `find_paths` - find all edges on paths between two nodes
- `exists_forward_path_between` / `exists_backward_path_between` - check path existence
- `path_from_root_to_node` / `path_from_node_to_node` - collect nodes along a path

## Utilities

- `assign_node_names` - assign auto-generated names (`NODE_0000000`, ...) to unnamed nodes
- `invert_edge` - reverse an edge's direction, updating adjacency lists on both endpoints
