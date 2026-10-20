# Reroot edge splitting accepts unbounded split fractions

The public reroot helper `fn split_edge()` accepts a raw `f64` split fraction. NaN and values outside $[0,1]$ give non-finite or negative branch lengths on the two new edges.

## Evidence

`fn split_edge()` [`packages/treetime-graph/src/reroot.rs#L50-L81`](../../packages/treetime-graph/src/reroot.rs#L50-L81) computes the parent-side length as `split_position * length` and the child-side length as `(1.0 - split_position) * length` at L63-L65, without a check of `split_position`. It also adds the new node at L56 before it reads the edge at L59. Every reroot calls it through `fn apply_root_at_edge()` [`packages/treetime/src/reroot/orchestrate.rs#L19`](../../packages/treetime/src/reroot/orchestrate.rs#L19) (call at L50).

On the reroot path, `RootTarget::on_edge()` [`packages/treetime/src/reroot/placement.rs#L29-L30`](../../packages/treetime/src/reroot/placement.rs#L29-L30) resolves the edge with `edge_endpoints()?` before `split_edge()` runs, so an unknown edge key gives an error before any mutation. The graph helpers resolve keys by index (`graph[edge_key]`, `graph[node_key]`), which panics on a missing key, as the project rule for node and edge data requires.

## Fix

- Add a finite `SplitFraction` type whose value $x$ satisfies $0\le x\le1$
- Use it in split evaluation, root-search results, and `split_edge()`
- Read the edge before the first mutation in `split_edge()`

## Validation

- Accept $0$, $0.5$, and $1$; reject negative, greater-than-one, NaN, and infinite inputs
- Split lengths are non-negative, finite, and sum to the original length
