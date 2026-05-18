# Pure transform architecture for graph traversals

Replace the current mutable-shared-state model (`Arc<RwLock<>>` on every node/edge payload, interior mutability during traversals) with a pure transform pattern where each traversal pass consumes typed input and produces typed output: `fn transform(input: &InputType) -> OutputType`.

The proposal is **not accepted** and **not implemented**. It records the motivation, the current architecture, the design dimensions, and the main tradeoffs.

## Motivation

The current graph architecture uses `Arc<RwLock<>>` on every payload. Any code holding a reference can read or write any field at any time. This is the Rust equivalent of v0's monkey-patched dynamic attributes: the type system enforces field existence but not access discipline.

Consequences:

- Locking overhead on every payload access, even in single-threaded sequential traversals
- No compile-time guarantee that a traversal pass only reads its inputs and writes its outputs
- Parallel BFS requires `RwLock` because callbacks mutate shared payloads in place
- Two independent lock acquisitions per node per parallel traversal step (one for graph payload, one for partition)
- Difficult to reason about data flow: any function receiving a `&Graph` can mutate any node/edge payload through `Arc::clone` + `write_arc()`
- `GraphNodeForward` and `GraphNodeBackward` reconstruct adjacency information from scratch per node per traversal step, acquiring O(degree) locks each time

A pure transform model would make data flow explicit at the type level: each pass declares what it reads and what it produces.

## Current architecture

### Graph payload plane

`struct Graph<N, E, D>` [packages/treetime-graph/src/graph.rs#L32](../../packages/treetime-graph/src/graph.rs#L32) stores nodes and edges behind `Arc<RwLock<>>`:

```rust
nodes: Vec<Option<Arc<RwLock<Node<N>>>>>,
edges: Vec<Option<Arc<RwLock<Edge<E>>>>>,
```

Each `struct Node<N>` [packages/treetime-graph/src/node.rs#L73](../../packages/treetime-graph/src/node.rs#L73) wraps its payload in a second `Arc<RwLock<N>>`:

```rust
pub struct Node<N: GraphNode> {
  key: GraphNodeKey,
  data: Arc<RwLock<N>>,       // payload: double-locked
  outbound_edges: Vec<GraphEdgeKey>,
  inbound_edges: Vec<GraphEdgeKey>,
  is_visited: AtomicBool,
}
```

Each `struct Edge<E>` [packages/treetime-graph/src/edge.rs#L71](../../packages/treetime-graph/src/edge.rs#L71) has the same pattern: `data: Arc<RwLock<E>>`.

Concrete payload types: `struct NodeAncestral` [packages/treetime/src/payload/ancestral.rs#L28](../../packages/treetime/src/payload/ancestral.rs#L28), `struct NodeTimetree` [packages/treetime/src/payload/timetree.rs#L17](../../packages/treetime/src/payload/timetree.rs#L17), `struct EdgeTimetree` [packages/treetime/src/payload/timetree.rs#L163](../../packages/treetime/src/payload/timetree.rs#L163). These carry tree-level data: names, time distributions, divergence, clock sets, branch length distributions, clock messages.

### Partition data plane

Partition structs (`struct PartitionMarginalSparse`, `struct PartitionMarginalDense`) store per-node and per-edge reconstruction state in `BTreeMap<GraphNodeKey, _>` and `BTreeMap<GraphEdgeKey, _>`, themselves behind `Arc<RwLock<Partition>>`. These carry profiles, messages, substitutions, indels, and likelihood contributions.

The two planes read from each other bidirectionally:

- Marginal passes read branch lengths from graph edge payloads, write profiles and messages to partition maps
- Branch-length optimization reads partition contributions, writes optimized lengths to graph edge payloads
- Timetree inference reads partition data to compute branch distributions, writes distributions to graph edge payloads
- Annotation reads partition substitutions, writes formatted mutation strings to graph node payloads

### Traversal callbacks

`struct GraphNodeForward` [packages/treetime-graph/src/graph_traverse.rs#L18](../../packages/treetime-graph/src/graph_traverse.rs#L18) and `struct GraphNodeBackward` [packages/treetime-graph/src/graph_traverse.rs#L96](../../packages/treetime-graph/src/graph_traverse.rs#L96) provide the callback context during traversals:

- `GraphNodeForward`: mutable access to current node payload + child edge payloads, read access to parent payloads
- `GraphNodeBackward`: mutable access to current node payload + parent edge payloads, read access to child payloads

Both types are reconstructed per node per traversal step by calling `graph.parents_of()` / `graph.children_of()`, each iterating edge keys, locking edges, reading source/target, and locking neighbor nodes. This is O(degree) lock acquisitions per node.

### Parallel BFS

The parallel BFS engine [packages/treetime-graph/src/breadth_first.rs#L139](../../packages/treetime-graph/src/breadth_first.rs#L139) processes frontiers via `rayon::into_par_iter()`. The explorer callback receives a `&SafeNodeRefMut<N>` (write guard on the node). The actual `GraphNodeForward`/`GraphNodeBackward` construction happens inside the callback wrapper at [packages/treetime-graph/src/graph_traverse.rs#L225-L253](../../packages/treetime-graph/src/graph_traverse.rs#L225-L253), acquiring additional locks on parents, children, and edges.

Partition locks are acquired inside the callback body. For example, `fn marginal_backward` [packages/treetime/src/commands/ancestral/marginal.rs#L85](../../packages/treetime/src/commands/ancestral/marginal.rs#L85) runs parallel BFS where each node callback acquires `partition.write_arc()` for every partition:

```rust
graph.par_iter_breadth_first_backward(|node| {
    // node already holds write_arc on Node payload
    for partition in partitions {
        let mut partition = partition.write_arc();  // second lock
        partition.process_node_backward(&node)?;
    }
});
```

### Topology mutations

Reroot (`fn split_edge` [packages/treetime/src/partition/algo/topology_cleanup/reroot.rs#L68](../../packages/treetime/src/partition/algo/topology_cleanup/reroot.rs#L68), `fn apply_reroot_topology` [packages/treetime/src/partition/algo/topology_cleanup/reroot.rs#L112](../../packages/treetime/src/partition/algo/topology_cleanup/reroot.rs#L112), `fn invert_edge` [packages/treetime-graph/src/edge.rs#L102](../../packages/treetime-graph/src/edge.rs#L102)), polytomy resolution (`fn resolve_polytomies_with_options` [packages/treetime/src/commands/timetree/optimization/polytomy.rs#L51](../../packages/treetime/src/commands/timetree/optimization/polytomy.rs#L51)), and edge collapse (`fn collapse_edge` [packages/treetime/src/partition/algo/topology_cleanup/collapse.rs#L34](../../packages/treetime/src/partition/algo/topology_cleanup/collapse.rs#L34)) mutate graph topology mid-pipeline via `graph.add_node()`, `graph.remove_edge()`, `graph.add_edge()`, `graph.collapse_edge()`. After topology changes, `fn reconcile_topology` patches partition maps with default entries for new nodes/edges.

### EM loop data flow

The timetree EM loop [packages/treetime/src/commands/timetree/run.rs#L244](../../packages/treetime/src/commands/timetree/run.rs#L244) cycles:

```
while !converged:
  1. run_refinement_iteration:
     a. update_marginal (partition backward + forward)
     b. optimize_branch_lengths (read partitions, write graph edges)
     c. reroot (topology mutation)
     d. run_timetree (read partitions + graph edges, write graph node/edge time data)
  2. record convergence metrics
```

Each step reads the previous step's mutations from the same mutable structures. The graph and partitions are both passed as `&mut`/`&` references with interior mutability allowing writes from anywhere.

## Design dimensions

### D1: Topology vs payload separation

Separate the immutable topology index (parent/child adjacency, node/edge keys, root/leaf sets) from mutable payload storage. Pre-compute a static adjacency structure once per topology epoch. Traversals consult the adjacency index directly instead of locking nodes to read edge keys.

Topology mutations (reroot, polytomy resolution, collapse) produce a new topology index. Payload storage is re-indexed to match.

### D2: Payload mutability model

Replace `Arc<RwLock<N>>` interior mutability with one of:

- Transform-based: traversal callback receives `&InputPayload`, returns `OutputPayload`. The traversal engine collects outputs into a new payload store. No mutation, no locks.
- Exclusive-ownership: traversal engine moves payloads out of the store, hands them to the callback as owned values, collects them back. Single owner at a time, no `Arc`, no `RwLock`, but allocation pattern unchanged.

Transform-based is cleaner for reasoning. Exclusive-ownership avoids the allocation of a separate output store.

### D3: Traversal callback signature

Current: `Fn(GraphNodeForward<N, E, D>) -> GraphTraversalContinuation` where `GraphNodeForward` contains `SafeNodePayloadRefMut<N>` (write guard) and `Vec<SafeEdgePayloadRefMut<E>>` (write guards on child edges).

Transform alternative:

```rust
Fn(TraversalContext<'_>) -> NodeOutput

struct TraversalContext<'a> {
    key: GraphNodeKey,
    is_root: bool,
    is_leaf: bool,
    node: &'a NodeInput,
    parent_edges: &'a [(GraphNodeKey, &'a EdgeInput)],
    child_edges: &'a [(GraphNodeKey, &'a EdgeOutput)],  // backward pass: children already computed
}
```

The callback reads from immutable references and returns a new value. The traversal engine places the return value into the output store.

### D4: Pipeline type staging

Current: same `NodeTimetree`/`EdgeTimetree` types throughout all passes. Fields are `Option<>` and populated at different stages. Which fields are valid depends on which passes have run.

Transform alternative: distinct types per pipeline stage, where each stage consumes the output of the previous stage:

```
FitchBackwardOutput -> FitchForwardOutput -> MarginalBackwardOutput -> MarginalForwardOutput
```

Each type contains exactly the fields that stage produces. Compile-time guarantee that downstream code cannot read fields that upstream has not yet computed.

### D5: Partition storage layout

Current: `BTreeMap<GraphNodeKey, _>` with O(log n) lookup per access.

Alternatives:

- Dense `Vec` indexed by `GraphNodeKey.as_usize()`. O(1) lookup, cache-friendly sequential access during traversals. Tombstones for removed nodes (matching the graph's existing `Vec<Option<>>` pattern).
- Flat arrays per field instead of per node. Structure-of-arrays layout for vectorizable traversal kernels.

### D6: Parallel BFS without locks

Current parallel BFS requires `RwLock` because callbacks mutate shared payloads in place. A transform model where each node produces output from immutable input is lock-free by construction: the traversal engine writes outputs to pre-allocated slots indexed by node key. Nodes at the same frontier level write to disjoint slots.

### D7: Dual data plane unification

Current architecture has two parallel mutable stores per node/edge (graph payload + partition map). A transform model could unify these into a single per-stage composite type, or keep them separate with explicit typed handoff between planes.

Unification simplifies reasoning (one store per node per stage). Separation preserves the ability to have multiple partitions per graph (multi-gene analysis) where each partition has its own per-node data but shares graph-level data.

A middle ground: graph-level data (topology, names, time distributions) as one typed layer, partition data as a separate typed layer, with explicit interface types for cross-plane reads.

### D8: Topology mutation strategy

Topology mutations (reroot, polytomy resolution, collapse) cannot produce "new output from old input" without copying the entire graph. Options:

- Epoch-based: topology mutations produce a new topology index plus a migration map (old key -> new key). Payload stores are re-indexed. Traversal passes within an epoch operate on immutable topology.
- Copy-on-mutate: topology mutations clone the graph, modify the clone, return it. Expensive for large trees but conceptually clean.
- In-place with fence: topology mutations are the one place that mutates in place, but a type-level fence prevents payload access during topology mutation (the graph is in a "topology-mutable" state that does not expose payload accessors).

Epoch-based is the natural fit: topology changes are infrequent (once per EM iteration), traversal passes are frequent (multiple per iteration).

## Existing partial precedents

The codebase already uses the transform pattern in limited scope:

- `Arc<Distribution>` values on `struct EdgeTimetree` are computed fresh each pass and stored as immutable shared references. `fn propagate_distributions_backward` [packages/treetime/src/commands/timetree/inference/backward_pass.rs#L17](../../packages/treetime/src/commands/timetree/inference/backward_pass.rs#L17) computes a new `Distribution`, wraps it in `Arc`, and stores it via `set_time_distribution()`. The distribution itself is never mutated after creation.
- `trait PartitionMarginalOps` [packages/treetime/src/partition/traits.rs#L102](../../packages/treetime/src/partition/traits.rs#L102) separates the graph (`&Graph`) from the partition (`&mut self`), so partition mutations do not require graph mutation. This is a step toward separating input (graph) from output (partition).
- `trait BranchTopology` [packages/treetime/src/partition/traits.rs#L26](../../packages/treetime/src/partition/traits.rs#L26) provides a read-only topology abstraction, decoupling partition operations from the concrete graph type.

## Impact assessment

### Scope

The change touches:

- `treetime-graph` crate: `struct Graph`, `struct Node`, `struct Edge`, all traversal methods, parallel BFS engine
- `treetime` crate: all partition types, all payload types, all traversal callbacks in all commands (ancestral, optimize, clock, timetree, mugration, prune)
- Test code: every test constructing a graph or using traversal callbacks

### Benefits

- Data flow visible in types: each pass's inputs and outputs are explicit
- No locking overhead for payload access
- Parallel BFS without `RwLock` contention
- Pre-computed adjacency index eliminates per-node O(degree) lock acquisitions during traversals
- Pipeline staging catches "read before write" bugs at compile time
- Dense `Vec` storage improves cache locality over `BTreeMap`

### Costs

- Large refactor touching every command and test
- Transform-based passes allocate a new output store per pass (mitigation: arena allocators or pre-allocated buffers)
- Pipeline type staging increases the number of types (mitigation: derive macros or type aliases for common combinations)
- Topology mutations require explicit re-indexing logic

### Risks

- Incremental migration is difficult because the current traversal API is pervasive
- Performance regression if output allocation dominates (unlikely for tree-sized data, needs measurement)
- Over-engineering if the locking overhead is not actually a bottleneck (needs profiling first)

## Validation plan

1. Profile the current architecture under large datasets (dengue/2000, sc2) to quantify locking overhead and `BTreeMap` lookup cost
2. Prototype D1 (topology/payload separation) and D5 (dense `Vec` storage) in isolation as they are independently valuable and measurable
3. Prototype D3 (transform callback signature) on a single traversal pass (marginal backward) to validate ergonomics
4. If prototypes show measurable improvement, plan incremental migration command by command

## Related documents

- [Graph-based phylogenetic representation](../decisions/graph-based-phylogenetic-representation.md) - documents the current `Graph<N, E, D>` design, `Arc<RwLock<>>` storage, typed payloads, parallel BFS, and the forward-looking rationale for DAG support
- [Partition system architecture](../decisions/partition-system-architecture.md) - documents the separation of tree topology from per-partition reconstruction state and the trait-based dispatch model
- [Dense and sparse sequence representation](../decisions/sequence-representation-dense-sparse.md) - documents the dense/sparse duality and trait-object interchangeability
- [Graph traversal algorithms](../algo/graph.md) - parallel BFS, DFS, path finding, edge contraction
- [Command modules contain shared operations that belong in domain layers](../issues/H-core-command-module-shared-ops-entanglement.md) - structural entanglement between command modules and shared operations (related: cleaner data flow would reduce cross-command coupling)
- [Dense and sparse partition types have structural and naming asymmetries](../issues/N-representation-dense-sparse-partition-asymmetry.md) - partition type asymmetries that a unified transform pipeline could resolve
- [Sparse marginal passes still use remove/insert pattern](../issues/N-ancestral-sparse-remove-insert-pattern.md) - in-place mutation pattern in sparse passes that a transform model would replace
- [Branch-length optimization: convergence and robustness](optimize-convergence-and-robustness.md) - documents the E-step/M-step loop where graph and partition data flow bidirectionally
