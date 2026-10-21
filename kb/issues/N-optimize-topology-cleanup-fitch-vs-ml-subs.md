# Topology cleanup in optimize loop uses Fitch substitutions instead of marginal MAP substitutions

The topology cleanup of each optimize iteration has two steps in `fn prune_and_merge_in_loop()`: the collapse of zero-optimal internal edges, then `fn resolve_polytomies()` (sibling merge, reversion hoist, and retirement of mutation-free helper nodes). Both steps decide from `fitch_subs()`, the Fitch substitutions of the sparse edge, not from the marginal MAP substitutions of the current ancestral reconstruction.

- `fn find_zero_optimal_internal_edges()` keeps an edge out of the collapse when any partition reports `edge_has_fitch_mutations()`, which reads `fitch_subs()` and the edge indels. Dense partitions always report no mutations
- `fn resolve_polytomies()` groups sibling edges by shared `fitch_subs()` (`fn merge_single_polytomy()`), counts reversions from `fitch_subs()` (`fn count_child_reversions()`), and retires helper nodes whose edges have empty `fitch_subs()`

The marginal MAP substitutions of a sparse partition are stored apart from the edge, in `SparseReconstruction::edges.estimates`, which the forward pass fills (`subs_ml`). At the point where `find_zero_optimal_internal_edges()` runs, the marginal update of the same iteration has filled them for every edge. `fitch_subs()` and the MAP substitutions can differ at ambiguous positions where marginal inference resolves differently than parsimony.

## Context

`fn collapse_edge()` composes the Fitch substitutions of the removed edge onto each child edge (`chain_fitch_subs()`). A composed edge has no matching marginal estimate until the next marginal update, so polytomy resolution after a collapse can use the MAP substitutions only for edges that the collapse did not change. The next iteration runs the marginal update again on the new topology.

## Suspected problems

- Zero-optimal selection on parsimony data: an edge where Fitch assigned a substitution at an ambiguous position, but the current MAP state has none, is kept out of the collapse. This is a conservative error (keeps edges that could be collapsed), not a correctness error (never collapses edges that carry current mutations)
- Polytomy resolution on mixed representations: for untouched sibling edges, `fitch_subs()` may differ from the MAP state at ambiguous sites. This can cause missed merges (shared in MAP, not in Fitch) or wrong merges (shared in Fitch, not in MAP), and different reversion counts for the hoist

## Investigation tasks

- Determine whether reading the MAP substitutions in `find_zero_optimal_internal_edges()` produces different zero-optimal edge sets on real datasets (compare the edge sets from `fitch_subs()` and from the marginal estimates after the marginal update)
- Determine whether reading the MAP substitutions for untouched sibling edges in `fn resolve_polytomies()` produces different merge and hoist decisions. Composed edges keep `fitch_subs()`, because they have no marginal estimate
- Measure the effect on iteration count and final topology

## Locations

- Zero-optimal edge identification: [`packages/treetime/src/optimize/run_loop.rs#L312-L331`](../../packages/treetime/src/optimize/run_loop.rs#L312-L331) (`find_zero_optimal_internal_edges`)
- Collapse and polytomy resolution: [`packages/treetime/src/optimize/run_loop.rs#L333-L385`](../../packages/treetime/src/optimize/run_loop.rs#L333-L385) (`prune_and_merge_in_loop`)
- Fitch mutation check: [`packages/treetime/src/partition/marginal/reconstruction.rs#L265-L274`](../../packages/treetime/src/partition/marginal/reconstruction.rs#L265-L274) (`edge_has_fitch_mutations`)
- Edge collapse: [`packages/treetime/src/optimize/topology/collapse.rs#L7-L43`](../../packages/treetime/src/optimize/topology/collapse.rs#L7-L43)
- Polytomy resolution: [`packages/treetime/src/optimize/topology/resolve_polytomy.rs#L17-L47`](../../packages/treetime/src/optimize/topology/resolve_polytomy.rs#L17-L47), with the sibling merge in [`packages/treetime/src/optimize/topology/merge_shared_mutations.rs`](../../packages/treetime/src/optimize/topology/merge_shared_mutations.rs) and the reversion count in [`packages/treetime/src/optimize/topology/hoist_reversions.rs#L15-L42`](../../packages/treetime/src/optimize/topology/hoist_reversions.rs#L15-L42)
- Fitch accessors: [`packages/treetime/src/partition/storage/sparse.rs#L105-L125`](../../packages/treetime/src/partition/storage/sparse.rs#L105-L125) (`SparseEdgeObs`)
- Marginal MAP substitutions: [`packages/treetime/src/partition/marginal/sparse/forward.rs#L67`](../../packages/treetime/src/partition/marginal/sparse/forward.rs#L67)
