# Graph traversal visit state is shared across queries

Breadth-first traversal stores session state on each persistent node through `AtomicBool is_visited` [`packages/treetime-graph/src/node.rs#L78`](../../packages/treetime-graph/src/node.rs#L78). Traversal marks nodes and public wrappers later reset the entire graph [`packages/treetime-graph/src/breadth_first.rs#L165-L191`](../../packages/treetime-graph/src/breadth_first.rs#L165-L191) [`packages/treetime-graph/src/graph_traverse.rs#L240-L268`](../../packages/treetime-graph/src/graph_traverse.rs#L240-L268).

`fn exists_path_between()` invokes the low-level traversal without resetting it, so a boolean query leaves hidden state behind [`packages/treetime-graph/src/find_paths.rs#L89`](../../packages/treetime-graph/src/find_paths.rs#L89). `find_paths()` compensates with explicit resets after each query [`packages/treetime-graph/src/find_paths.rs#L29-L39`](../../packages/treetime-graph/src/find_paths.rs#L29-L39).

## State ownership

`Node<N>` combines persistent graph identity, payload, inbound and outbound edges, and the transient `is_visited` flag. Public traversal methods accept `&self`, so their signatures appear read-only even though they mutate every visited node and later perform a graph-wide reset.

The marker has no traversal identifier. It can represent only `visited by some traversal`, which is insufficient when queries overlap or when a low-level caller exits without the expected reset.

## Failure mode

- Sequential query results depend on whether callers remembered the reset protocol.
- Concurrent traversals can observe each other's visited nodes.
- One traversal can clear markers while another traversal still needs them.

Atomic access prevents a memory race on one bit; it does not distinguish traversal sessions. Visit state must be local to a traversal, or exclusivity must be represented and enforced by the API.

## Required design

- O1. Store visit state in the traversal operation, keyed by `GraphNodeKey`. This makes independent traversals composable and keeps nodes limited to persistent graph state.
- O2. Enforce graph-wide traversal exclusivity with an explicit guard and make every traversal path use it. This preserves shared markers but serializes traversal and exposes the exclusivity contract.

O1 matches the current `&Graph` query surface and avoids hidden mutation. Any implementation must preserve deterministic traversal order and cycle handling.

## Validation

- Sequential path queries require no reset and remain independent.
- Two simultaneous forward/backward traversals of one graph produce the same results as isolated traversals.
- Early errors and cancellation leave no graph state that affects a later query.
- Chain, multifurcation, forest, and cyclic fixtures exercise visit semantics.

## Related issues

- [H-graph-lock-topology-is-exposed-to-consumers.md](H-graph-lock-topology-is-exposed-to-consumers.md)
