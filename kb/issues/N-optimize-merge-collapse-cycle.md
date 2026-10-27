# Optimize loop alternates between merging and collapsing the same siblings

The topology cleanup of each optimize iteration has two steps in `fn prune_and_merge_in_loop()`: the collapse of zero-optimal internal edges, then `fn resolve_polytomies()` (sibling merge, reversion hoist, and retirement of mutation-free helper nodes). The merge groups siblings by Fitch substitutions that the marginal estimates of the iteration also place on each of them. Some of these merges create an edge whose ML length is zero, so the next iteration collapses it and the one after merges the siblings again.

## Context

- The collapse does not read substitutions: `fn find_zero_optimal_internal_edges()` collapses every internal edge whose ML length is zero, unless an indel term of the likelihood requires a positive length
- The sibling merge (`fn build_mutation_index()`) and the reversion count of the hoist (`fn count_child_reversions()`) read the Fitch substitutions that the marginal estimate of the edge also contains (`fn SparseReconstruction::supported_fitch_subs()`). An edge without an estimate, as in `prune`, keeps all its Fitch substitutions
- The rewrites keep the estimates of the edges they change: the collapse composes them, the merge and the hoist split them like the Fitch substitutions. The root slide of the hoist (`fn slide_bifurcating_root_for_child()`) drops the estimates of the two root edges
- `fn run_optimize_loop()` does not stop on convergence or on a decrease in the iteration after a topology change, because the new edges still have their initial lengths; it stops on an oscillation at once. The best iteration keeps its tree, names, and partitions, and the loop returns the tree with the highest likelihood

On flu/h3n2/500, a merge of 7 siblings that share `T906C` in both the Fitch and the marginal substitutions creates an edge that the next iteration optimizes to zero and collapses. The loop stops on the oscillation and returns the tree without the merge, which has the 205 splits of v0. On sc2/5100, iterations 4 to 7 alternate between collapsing 3 edges and merging 3 sibling groups, and the loop returns iteration 5.

## Suspected problems

- The marginal estimates are the most likely state of each node on its own. A substitution can be the most likely on each of several sibling edges and still not on the edge that groups them, so the merge proposes a change that the likelihood rejects. The loop ends on the best tree, but spends 2 iterations per such cycle
- `fn retire_created_helpers()` still retires helper nodes by empty `fitch_subs()`

## Investigation tasks

- Decide whether a merge should be accepted only when it does not lower the likelihood, for example by scoring the grouped edge before the merge
- Decide whether helper retirement should read the marginal estimates

## Locations

- Loop and stop rule: [`packages/treetime/src/optimize/run_loop.rs`](../../packages/treetime/src/optimize/run_loop.rs) (`run_optimize_loop`, `stop_reason`, `BestIteration`)
- Edge collapse: [`packages/treetime/src/optimize/topology/collapse.rs`](../../packages/treetime/src/optimize/topology/collapse.rs)
- Polytomy resolution: [`packages/treetime/src/optimize/topology/resolve_polytomy.rs`](../../packages/treetime/src/optimize/topology/resolve_polytomy.rs), with the sibling merge in [`packages/treetime/src/optimize/topology/merge_shared_mutations.rs`](../../packages/treetime/src/optimize/topology/merge_shared_mutations.rs) and the hoist in [`packages/treetime/src/optimize/topology/hoist_reversions.rs`](../../packages/treetime/src/optimize/topology/hoist_reversions.rs)
- Supported substitutions: [`packages/treetime/src/partition/marginal/reconstruction.rs`](../../packages/treetime/src/partition/marginal/reconstruction.rs) (`supported_fitch_subs`)
