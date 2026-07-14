# Reroot edge splitting lacks validation and failure atomicity

Public reroot helpers accept raw split fractions and mutate the graph before validating all keys. Type-valid inputs can create non-finite or negative branch lengths, and stale keys can panic after leaving an orphan node.

## Evidence

`fn split_edge()` [`packages/treetime-graph/src/reroot.rs#L71`](../../packages/treetime-graph/src/reroot.rs#L71) inserts a node before resolving the edge and uses `expect`. Its `f64` split accepts NaN and values outside $[0,1]$.

## Potential solutions

- O1. Validate a bounded split type and all keys, then commit a completely constructed topology delta.
- O2. Mutate incrementally with rollback guards. This requires every future fallible mutation to participate in rollback correctly.

## Recommendation

Parse a `SplitFraction` type whose value $x$ is finite and satisfies $0\le x\le1$. Resolve all node and edge keys before mutation, construct the topology change completely, and commit only after every fallible step succeeds.

## Related issues

- [N-reroot-tip-resolution-untested-errors.md](N-reroot-tip-resolution-untested-errors.md)
