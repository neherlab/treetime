# Reversion hoist in optimize loop uses Fitch substitutions instead of marginal MAP substitutions

The topology cleanup of each optimize iteration has two steps in `fn prune_and_merge_in_loop()`: the collapse of zero-optimal internal edges, then `fn resolve_polytomies()` (sibling merge, reversion hoist, and retirement of mutation-free helper nodes). The reversion hoist and the helper retirement decide from `fitch_subs()`, the Fitch substitutions of the sparse edge, not from the marginal MAP substitutions of the current ancestral reconstruction.

- `fn count_child_reversions()` counts reversions from `fitch_subs()`, and `fn retire_created_helpers()` retires helper nodes whose edges have empty `fitch_subs()`

The marginal MAP substitutions of a sparse partition are stored apart from the edge, in `SparseReconstruction::edges.estimates`, which the forward pass fills (`subs_ml`). The marginal update of each iteration fills them for every edge. `fitch_subs()` and the MAP substitutions can differ at ambiguous positions where marginal inference resolves differently than parsimony.

## Context

The other two steps no longer act on a Fitch placement that the marginal reconstruction does not share:

- The collapse does not read substitutions: `fn find_zero_optimal_internal_edges()` collapses every internal edge whose ML length is zero, unless an indel term of the likelihood requires a positive length. The collapse composes the Fitch substitutions, the indels, and the marginal estimates of the removed edge onto the child edges; a child edge gets no estimate when either edge lacks one
- The sibling merge (`fn build_mutation_index()`) groups siblings by the Fitch substitutions that the marginal estimate of the edge also contains (`fn SparseReconstruction::supports_fitch_sub()`); an edge without an estimate, as in `prune` or on an edge a hoist created, keeps all its Fitch substitutions. The merge splits the estimates of the grouped edges like their Fitch substitutions

On flu/h3n2/500, the Fitch veto on the collapse left 18 of 223 internal edges at length zero. Without the veto but with a Fitch merge, the loop collapsed such an edge, moved its Fitch substitution onto the children, merged the children again, and stopped with 12 extra internal edges, most at the initial length of one mutation and without a marginal substitution. With the marginal merge, 6 extra internal edges remain against v0 (211 against 205), and 1 with `--no-hoist-reversions`.

## Suspected problems

- Hoist on parsimony data: a reversion that Fitch places but the marginal reconstruction does not can trigger a hoist. With the hoist, 5 more internal edges remain on flu/h3n2/500, 4 of them original edges at about half a mutation without a marginal substitution
- A merge taken from the estimates of one iteration can lower the likelihood of the next. On flu/h3n2/500 without the hoist, the 1 remaining extra edge comes from such a merge (`T906C` shared by 7 siblings); the loop then stopped because the likelihood decreased, so the new edge kept its initial length

## Investigation tasks

- Read the MAP substitutions in `fn count_child_reversions()` and keep the hoist data consistent the way the merge does
- Decide whether the loop should re-optimize edges created by the last topology change before it stops

## Locations

- Collapse and polytomy resolution: [`packages/treetime/src/optimize/run_loop.rs`](../../packages/treetime/src/optimize/run_loop.rs) (`prune_and_merge_in_loop`)
- Edge collapse: [`packages/treetime/src/optimize/topology/collapse.rs`](../../packages/treetime/src/optimize/topology/collapse.rs)
- Polytomy resolution: [`packages/treetime/src/optimize/topology/resolve_polytomy.rs`](../../packages/treetime/src/optimize/topology/resolve_polytomy.rs), with the sibling merge in [`packages/treetime/src/optimize/topology/merge_shared_mutations.rs`](../../packages/treetime/src/optimize/topology/merge_shared_mutations.rs) and the reversion count in [`packages/treetime/src/optimize/topology/hoist_reversions.rs`](../../packages/treetime/src/optimize/topology/hoist_reversions.rs)
- Fitch accessors: [`packages/treetime/src/partition/storage/sparse.rs`](../../packages/treetime/src/partition/storage/sparse.rs) (`SparseEdgeObs`)
- Marginal MAP substitutions: [`packages/treetime/src/partition/marginal/sparse/forward.rs`](../../packages/treetime/src/partition/marginal/sparse/forward.rs)
