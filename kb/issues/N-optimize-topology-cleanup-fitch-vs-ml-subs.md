# Polytomy resolution in optimize loop uses Fitch substitutions instead of marginal MAP substitutions

The topology cleanup of each optimize iteration has two steps in `fn prune_and_merge_in_loop()`: the collapse of zero-optimal internal edges, then `fn resolve_polytomies()` (sibling merge, reversion hoist, and retirement of mutation-free helper nodes). The second step decides from `fitch_subs()`, the Fitch substitutions of the sparse edge, not from the marginal MAP substitutions of the current ancestral reconstruction.

- `fn resolve_polytomies()` groups sibling edges by shared `fitch_subs()` (`fn merge_single_polytomy()`), counts reversions from `fitch_subs()` (`fn count_child_reversions()`), and retires helper nodes whose edges have empty `fitch_subs()`

The marginal MAP substitutions of a sparse partition are stored apart from the edge, in `SparseReconstruction::edges.estimates`, which the forward pass fills (`subs_ml`). The marginal update of each iteration fills them for every edge. `fitch_subs()` and the MAP substitutions can differ at ambiguous positions where marginal inference resolves differently than parsimony.

## Context

The collapse step no longer reads substitutions: `fn find_zero_optimal_internal_edges()` collapses every internal edge whose ML length is zero, unless an indel term of the likelihood requires a positive length. With Fitch substitutions as a veto, an edge whose ambiguous change Fitch had placed on it stayed in the tree at length zero; on flu/h3n2/500 these were 18 of 223 internal edges. The collapse composes the substitutions and indels of the removed edge onto the child edges, so nothing is lost.


`fn collapse_edge()` composes the Fitch substitutions of the removed edge onto each child edge (`chain_fitch_subs()`). A composed edge has no matching marginal estimate until the next marginal update, so polytomy resolution after a collapse can use the MAP substitutions only for edges that the collapse did not change. The next iteration runs the marginal update again on the new topology.

## Suspected problems

- Polytomy resolution on mixed representations: for untouched sibling edges, `fitch_subs()` may differ from the MAP state at ambiguous sites. This can cause missed merges (shared in MAP, not in Fitch) or wrong merges (shared in Fitch, not in MAP), and different reversion counts for the hoist

## Investigation tasks

- Determine whether reading the MAP substitutions for untouched sibling edges in `fn resolve_polytomies()` produces different merge and hoist decisions. Composed edges keep `fitch_subs()`, because they have no marginal estimate
- Measure the effect on iteration count and final topology

## Locations

- Collapse and polytomy resolution: [`packages/treetime/src/optimize/run_loop.rs#L333-L385`](../../packages/treetime/src/optimize/run_loop.rs#L333-L385) (`prune_and_merge_in_loop`)
- Edge collapse: [`packages/treetime/src/optimize/topology/collapse.rs#L7-L43`](../../packages/treetime/src/optimize/topology/collapse.rs#L7-L43)
- Polytomy resolution: [`packages/treetime/src/optimize/topology/resolve_polytomy.rs#L17-L47`](../../packages/treetime/src/optimize/topology/resolve_polytomy.rs#L17-L47), with the sibling merge in [`packages/treetime/src/optimize/topology/merge_shared_mutations.rs`](../../packages/treetime/src/optimize/topology/merge_shared_mutations.rs) and the reversion count in [`packages/treetime/src/optimize/topology/hoist_reversions.rs#L15-L42`](../../packages/treetime/src/optimize/topology/hoist_reversions.rs#L15-L42)
- Fitch accessors: [`packages/treetime/src/partition/storage/sparse.rs#L105-L125`](../../packages/treetime/src/partition/storage/sparse.rs#L105-L125) (`SparseEdgeObs`)
- Marginal MAP substitutions: [`packages/treetime/src/partition/marginal/sparse/forward.rs#L67`](../../packages/treetime/src/partition/marginal/sparse/forward.rs#L67)
