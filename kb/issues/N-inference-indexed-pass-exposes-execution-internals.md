# Indexed inference passes expose scheduling and storage internals

`IndexedPass` reuses graph-generated forward and backward frontiers, but its interface exposes range offsets, completed/future slices, indexed payload storage, and publication behavior to every inference callback [`packages/treetime/src/partition/indexed_pass.rs#L75`](../../packages/treetime/src/partition/indexed_pass.rs#L75) [`packages/treetime/src/partition/indexed_pass.rs#L175`](../../packages/treetime/src/partition/indexed_pass.rs#L175).

## Impact

Fitch, marginal, clock, and timetree algorithms must understand execution ordering and indexed storage mechanics in addition to their domain operation. Failure and cancellation semantics therefore remain coupled to callback implementations even though topology planning is shared.

`with_indexed_graph_payloads()` removes payload maps from the graph, invokes a callback, restores the maps, and then propagates the callback error. Mutations made before an error are therefore restored as partially committed state rather than rolled back. This failure contract is tracked separately as a high-severity graph defect.

## Design question

Determine whether the execution module should own dependency-safe node access, cancellation, and commit/rollback while domain code supplies only node work, or whether the low-level interface is required for a measured performance property.

## Required execution contract

- Own topology planning, indexed storage, forward/backward execution, cancellation, and publication in one module.
- Present domain callbacks with the completed dependencies and writable current output needed for one node.
- Define whether an error rolls back all payload changes or commits a documented complete prefix.
- Restore graph payload ownership on every success, error, and cancellation path.

## Validation

- Chain, balanced tree, multifurcation, forest, and malformed-topology fixtures exercise dependency ordering.
- One-worker and multi-worker execution produce equivalent published payloads.
- Injected callback errors at each frontier verify the selected atomicity contract.
- Fitch, marginal, clock, and timetree tests assert domain outcomes without depending on storage ranges.

## Related issues

- [H-graph-indexed-payload-extraction-exposes-invalid-state.md](H-graph-indexed-payload-extraction-exposes-invalid-state.md)
- [M-inference-fallible-parallel-passes-partially-commit.md](M-inference-fallible-parallel-passes-partially-commit.md)
