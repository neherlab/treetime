# Make indexed payload passes failure-atomic with scratch state

Prevent shared graph readers from observing default payloads and commit indexed-pass results only after successful completion.

## Required changes

- Apply the scratch-state design to every pass using `with_indexed_graph_payloads()`. These passes read stable graph topology and transform node and edge payloads; they do not intrinsically mutate topology.
- Snapshot the complete valid node and edge payload maps without replacing published graph payloads with defaults.
- Run the indexed pass entirely against independent scratch maps. On callback error or panic, discard the scratch maps and leave every published payload unchanged.
- Add a graph-level payload transaction gate. Every logical operation that reads or writes one or more published payloads must hold its typed shared guard for the operation's complete lifetime; direct unguarded multi-payload access is no longer part of the API. Audit all payload callers and migrate them to the guard.
- After successful completion and structural validation, acquire the transaction gate exclusively, acquire node and edge payload locks in deterministic key order, recheck the snapshot generation, and replace every payload before releasing the exclusive guard. A stale generation or lock/precondition failure is an internal error and must occur before the first replacement, leaving the published graph unchanged.
- Remove the extraction/restoration path from this helper; an unwind guard is neither an alternative nor part of this ticket's implementation.
- Document the invariant that readers observe either the complete pre-pass payload state or the complete post-pass payload state, never defaults or a partial mixture.

## Validation

- Concurrent reader/writer tests using the public guarded API, proving readers observe only a complete pre-pass or post-pass state and writers cannot deadlock when locking payloads in key order.
- Inject worker error, panic, structural validation failure, and commit precondition failure; compare the complete graph with its pre-pass state.
- Deterministic successful result under one and multiple workers.
- Full lint and test suite.

## Related issues

- Source: [kb/issues/H-graph-indexed-payload-extraction-exposes-invalid-state.md](../issues/H-graph-indexed-payload-extraction-exposes-invalid-state.md)
