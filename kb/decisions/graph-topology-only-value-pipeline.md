# Graph is topology only; inference data flows as value maps

TreeTime v1 stores phylogenetic structure and inference results separately. The directed graph (`struct Graph` in the `treetime-graph` crate) carries topology only: a node holds an identity key and its inbound and outbound edge keys, an edge holds an identity key and its source and target node keys. The graph has no node payload, no edge payload, and no graph-level data slot, and it is not generic over stored data.

Every value an inference stage reads or writes lives outside the graph, addressed by the graph's stable keys:

- Per-node and per-edge data flow as maps keyed by `GraphNodeKey` and `GraphEdgeKey`: node names, branch lengths, ancestral sequences and profiles, mutation lists, date and time distributions, node divergence, and outlier flags.
- Formerly graph-global data flows as explicit pass arguments and per-command result values: the substitution (GTR) model, the per-gene partitions, the clock model, and the clock-regression results. The per-command result types own them (for example `AncestralOutput`, `OptimizeOutput`, `PruneOutput`, `ClockOutput`, `TimetreeOutput`).

Each inference stage is a function from declared inputs to declared outputs. A stage receives the value maps it reads and returns new maps; it does not reach through a graph reference to read or mutate shared state.

## Context

An earlier design made the graph generic over a node payload, an edge payload, and a graph-level data type (`Graph<N, E, D>`), and stored each payload behind an `Arc<RwLock<...>>`. The per-gene partition was a second shared, lock-guarded store. Lock guards protected concurrent access, but they did not restrict a stage to its declared inputs and outputs: any stage holding a graph reference could read or write any payload field, and which fields a stage required or produced was implicit.

## Decision

- The graph carries topology only: no node, edge, or graph-level payload, and no data type parameter.
- Per-node and per-edge inference data is held in value maps keyed by the node and edge keys, owned by the pipeline that runs the stages.
- Graph-level inputs and results (GTR model, partitions, clock model, regression results) pass as arguments and return as values held by the per-command result types.
- Inference stages are value transforms. Parallel passes operate on owned per-node and per-edge values and take no payload lock.

This amends [graph-based-phylogenetic-representation.md](graph-based-phylogenetic-representation.md): the directed-graph model and its tree-and-network rationale stand, but the graph no longer carries typed node, edge, or graph payloads. Branch-specific data keeps the edge as its owner by being keyed to the edge, now in a side table rather than in an edge field.

## Alternatives considered

- Keep a typed graph-level data slot (`Graph<D>`) for graph-global values such as the GTR model and partitions. Rejected: a shared mutable field on the graph is the same implicit channel this decision removes. A stage can read or write it outside its declared contract, and it prevents lock-free parallel passes.
- Keep typed node and edge payloads behind locks. Rejected for the same reason: access runs through a shared graph reference rather than through declared stage inputs and outputs.

## Consequences

- `--output-tree-graph-json` serializes the graph, which is now topology only: nodes, edges, roots, and leaves. Values that the previous graph dump embedded through its data slot (for example the substitution-model name and the alignment mask) are no longer in the graph-JSON file. They remain available through their dedicated outputs (the GTR model file and the node-data and tree outputs). This is a deliberate change to the content of that one output.
- All other command outputs (ancestral, clock, mugration, optimize, prune, timetree) are byte-identical to the previous tree-and-payload representation, confirmed across the output-comparison matrix and a thread-count determinism sweep.
- A stage's inputs and outputs are explicit in its signature, and the type of each value map names the stage that produces it.
- Parallel passes need no payload lock. The graph's node and edge handles remain reference-counted, read-locked structures for topology, written only during single-threaded graph construction and rerooting.

## Related decisions

- [graph-based-phylogenetic-representation.md](graph-based-phylogenetic-representation.md): the directed-graph topology model.
- [partition-system-architecture.md](partition-system-architecture.md): separation of topology from per-partition state.
