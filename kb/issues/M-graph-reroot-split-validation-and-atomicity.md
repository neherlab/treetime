# Reroot edge splitting lacks validation and failure atomicity

Public reroot helpers accept raw split fractions and mutate the graph before validating all keys. Type-valid inputs can create non-finite or negative branch lengths, and stale keys can panic after leaving an orphan node.

## Evidence

`fn split_edge()` [`packages/treetime-graph/src/reroot.rs#L55`](../../packages/treetime-graph/src/reroot.rs#L55) inserts a node at L61 before resolving the edge, and uses `expect` at L64. Its `f64` split accepts NaN and values outside $[0,1]$. The clock reroot [`packages/treetime/src/clock/reroot.rs#L146`](../../packages/treetime/src/clock/reroot.rs#L146) and the generic reroot orchestration [`packages/treetime/src/reroot/orchestrate.rs#L118`](../../packages/treetime/src/reroot/orchestrate.rs#L118) call it.

`fn remove_node_if_trivial()` [`packages/treetime-graph/src/reroot.rs#L111`](../../packages/treetime-graph/src/reroot.rs#L111) and `fn trivial_node_branch_lengths()` [`packages/treetime-graph/src/reroot.rs#L152`](../../packages/treetime-graph/src/reroot.rs#L152) also resolve keys with `expect`.

## Potential solutions

- O1. Validate a bounded split type and all keys, then commit a completely constructed topology delta.
- O2. Mutate incrementally with rollback guards. This requires every future fallible mutation to participate in rollback correctly.

## Recommendation

Parse a `SplitFraction` type whose value $x$ is finite and satisfies $0\le x\le1$. Resolve all node and edge keys before mutation, construct the topology change completely, and commit only after every fallible step succeeds.

## Fix (O1)

- Add a finite `SplitFraction` type bounded to $[0,1]$
- Use it in split evaluation, root-search results, and the graph mutation APIs
- Resolve every supplied key and every required edge and node before the first mutation
- Replace recoverable `expect` calls with contextual project errors
- Commit the complete edge split atomically

## Validation

- Accept $0$, $0.5$, and $1$; reject negative, greater-than-one, NaN, and infinite inputs
- Split lengths are non-negative, finite, and sum to the original length
- Inject stale node and edge keys and compare the entire graph before and after the error
