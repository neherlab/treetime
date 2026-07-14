# Consolidate dependency-frontier schedulers

Replace the graph and partition frontier implementations with one generic scheduler.

## Required changes

- Model dependency counts and ready nodes once.
- Provide domain callbacks for work and deterministic publication.
- Centralize cancellation, worker error propagation, and incomplete-graph detection.
- Preserve existing forward/backward traversal order contracts.

## Validation

- Chain, balanced tree, multifurcation, forest, and cycle/error fixtures.
- One-worker and multi-worker equivalence.
- Inject worker errors and verify no further commits occur.
- Full lint and test suite.

## Related issues

- Source: [kb/issues/N-graph-dependency-frontier-schedulers-duplicated.md](../issues/N-graph-dependency-frontier-schedulers-duplicated.md)
