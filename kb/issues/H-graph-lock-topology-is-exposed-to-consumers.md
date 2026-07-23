# Graph lock topology is exposed to consumers

The graph API returns `Arc<RwLock<Node<N>>>` and `Arc<RwLock<Edge<E>>>`, while each node and edge returns another locked payload [`packages/treetime-graph/src/graph.rs#L13`](../../packages/treetime-graph/src/graph.rs#L13) [`packages/treetime-graph/src/node.rs#L71`](../../packages/treetime-graph/src/node.rs#L71) [`packages/treetime-graph/src/edge.rs#L67`](../../packages/treetime-graph/src/edge.rs#L67).

Ordinary domain operations therefore contain chains such as `edge.read_arc().payload().read_arc()`. These chains occur across I/O, rerooting, topology, clock, pruning, marginal inference, optimization, and timetree code. Consumers depend on wrapper nesting, lock granularity, graph element representation, and payload ownership.

## Blast radius

Structural inspection found the exact nested read chain across numerous production modules, with additional write and endpoint variants. Representative consumers include Auspice and Graphviz I/O, root-cost evaluation, marginal inference, timetree likelihood, and command output. A change to payload ownership or lock placement therefore requires coordinated edits across domain and format crates.

The chains also make lock scope difficult to review. A caller controls when each guard is created and dropped, and multi-object operations can acquire node, edge, and payload locks in inconsistent orders.

## Design question

Determine which graph-owned query and command operations can hide intermediate guards while preserving concurrency and performance requirements. Raw wrapper traversal should remain available only where the caller genuinely owns a multi-object lock protocol.

## Design axes

- Payload access: graph-owned read/write closures versus purpose-specific payload guards.
- Endpoint access: one operation that acquires an edge and its endpoint payloads versus separate caller-managed guards.
- Multi-object mutation: explicit transaction/ordering API versus exposed raw locks for a documented narrow set of algorithms.

The abstraction must avoid returning references beyond guard lifetimes and must preserve the parallel read/write behavior required by inference.

## Validation

- Lock-order tests and injected contention cover node/edge/payload combinations used by rerooting and inference.
- Domain operations no longer depend on `Node`, `Edge`, and payload wrapper nesting.
- Throughput benchmarks compare traversal and marginal passes before raw access is restricted.
- Poisoning/error paths retain original context and cannot leave partially applied graph commands.

## Related issues

- [H-graph-traversal-visit-state-is-shared.md](H-graph-traversal-visit-state-is-shared.md)
- [H-graph-indexed-payload-extraction-exposes-invalid-state.md](H-graph-indexed-payload-extraction-exposes-invalid-state.md)
