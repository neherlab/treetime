# Optimize loop keeps an unoptimized topology change when it stops

The topology cleanup of each optimize iteration has two steps in `fn prune_and_merge_in_loop()`: the collapse of zero-optimal internal edges, then `fn resolve_polytomies()` (sibling merge, reversion hoist, and retirement of mutation-free helper nodes). The merge and the hoist decide from the marginal estimates of the iteration. A topology change made from these estimates can lower the likelihood of the next iteration, and the loop then stops with the changed topology and the starting lengths of the new edges.

## Context

The cleanup steps no longer act on a Fitch placement that the marginal reconstruction does not share:

- The collapse does not read substitutions: `fn find_zero_optimal_internal_edges()` collapses every internal edge whose ML length is zero, unless an indel term of the likelihood requires a positive length
- The sibling merge (`fn build_mutation_index()`) and the reversion count of the hoist (`fn count_child_reversions()`) read the Fitch substitutions that the marginal estimate of the edge also contains (`fn SparseReconstruction::supported_fitch_subs()`). An edge without an estimate, as in `prune`, keeps all its Fitch substitutions
- The rewrites keep the estimates of the edges they change: the collapse composes them, the merge and the hoist split them like the Fitch substitutions. The root slide of the hoist (`fn slide_bifurcating_root_for_child()`) drops the estimates of the two root edges, which then fall back to their Fitch substitutions
- The sparse likelihood is computed from the Fitch substitutions, so the rewrites still move them; the marginal estimates only select which rewrite to make

On flu/h3n2/500, with Fitch decisions the loop left 18 internal edges at length zero (Fitch veto on the collapse), or collapsed such an edge and merged its children again, ending with 12 internal edges more than v0. With marginal decisions, the tree has the 205 splits of v0 and 1 more (206 internal edges), with the same 1507 marginal substitutions.

## Suspected problems

- The 1 remaining extra edge on flu/h3n2/500 comes from a merge of 7 siblings that share `T906C` in both the Fitch and the marginal substitutions. The likelihood of the next iteration is lower (-13101.01 to -13101.71), the loop stops on `ConvergenceReason::Worsened`, and the new edge keeps its starting length of one mutation. The reconstruction of the final tree places `C906T` on 6 lineages and no substitution on the new edge
- `fn retire_created_helpers()` still retires helper nodes by empty `fitch_subs()`

## Investigation tasks

- Decide whether the loop should re-optimize the edges of the last topology change before it stops, or revert a topology change after which the likelihood decreased
- Decide whether helper retirement should read the marginal estimates

## Locations

- Collapse and polytomy resolution: [`packages/treetime/src/optimize/run_loop.rs`](../../packages/treetime/src/optimize/run_loop.rs) (`prune_and_merge_in_loop`, `stop_reason`)
- Edge collapse: [`packages/treetime/src/optimize/topology/collapse.rs`](../../packages/treetime/src/optimize/topology/collapse.rs)
- Polytomy resolution: [`packages/treetime/src/optimize/topology/resolve_polytomy.rs`](../../packages/treetime/src/optimize/topology/resolve_polytomy.rs), with the sibling merge in [`packages/treetime/src/optimize/topology/merge_shared_mutations.rs`](../../packages/treetime/src/optimize/topology/merge_shared_mutations.rs) and the hoist in [`packages/treetime/src/optimize/topology/hoist_reversions.rs`](../../packages/treetime/src/optimize/topology/hoist_reversions.rs)
- Supported substitutions: [`packages/treetime/src/partition/marginal/reconstruction.rs`](../../packages/treetime/src/partition/marginal/reconstruction.rs) (`supported_fitch_subs`)
- Marginal MAP substitutions: [`packages/treetime/src/partition/marginal/sparse/forward.rs`](../../packages/treetime/src/partition/marginal/sparse/forward.rs)
