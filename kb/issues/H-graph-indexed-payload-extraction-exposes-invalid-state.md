# Indexed payload extraction exposes invalid shared graph state

## Problem

`with_indexed_graph_payloads()` removes payload maps from the graph with `mem::take`, releases the locks, and runs fallible work while the shared graph contains default payloads. Concurrent readers can observe invalid data; panic or restoration failure can leave the defaults permanently installed.

The extraction and restoration transaction is implemented in `fn with_indexed_graph_payloads()` [packages/treetime/src/partition/indexed_pass.rs#L8-L60](../../packages/treetime/src/partition/indexed_pass.rs#L8-L60). Marginal and timetree passes invoke fallible parallel frontier callbacks while those payloads are detached [packages/treetime/src/partition/marginal_passes.rs#L37-L88](../../packages/treetime/src/partition/marginal_passes.rs#L37-L88).

## Options

- O1. Hold exclusive ownership for the complete extraction/work/restoration transaction with unwind-safe restoration.
- O2. Compute against independent scratch payloads and atomically commit completed maps.

Recommendation: O2 when parallel work does not require graph mutation; it keeps the published graph valid throughout and gives failure atomicity. Use O1 only where exclusive graph mutation is intrinsic.

## Recommendation

Use O2 for `with_indexed_graph_payloads()` and every current caller: snapshot valid payloads, compute complete replacement maps in scratch state, then atomically commit only after successful computation and structural validation. The affected passes read stable topology, so no current caller needs the O1 extraction/restoration design. O1 remains an option only for a future operation that intrinsically mutates topology and must not be mixed into this ticket.

## Related issues

- [M-inference-fallible-parallel-passes-partially-commit.md](M-inference-fallible-parallel-passes-partially-commit.md)
- [N-graph-dependency-frontier-schedulers-duplicated.md](N-graph-dependency-frontier-schedulers-duplicated.md)
