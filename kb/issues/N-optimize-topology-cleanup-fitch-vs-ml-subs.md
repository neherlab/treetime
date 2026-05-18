# Topology cleanup in optimize loop uses Fitch substitutions instead of marginal MAP substitutions

The optimize loop's topology cleanup step (`find_zero_optimal_internal_edges` + `merge_shared_mutation_branches`) uses `fitch_subs()` to decide which edges carry mutations and which sibling edges share mutations. `fitch_subs()` contains the compression-time Fitch-algorithm substitutions, not the post-marginal MAP-derived substitutions that reflect the current ancestral reconstruction.

At the point where `find_zero_optimal_internal_edges` runs, `update_marginal()` has already populated `subs_ml` with the current MAP state. For non-collapsed edges, `fitch_subs()` and `subs_ml` can differ at ambiguous positions where marginal inference resolves differently than parsimony. The topology cleanup decisions are therefore based on the parsimony-era representation, not the current maximum a posteriori state.

## Context

This is the documented current behavior, not a regression. Git archaeology shows that `find_zero_optimal_internal_edges` always read the raw `subs` field (now `subs_fitch` behind the `fitch_subs()` accessor). The `edge_subs()` trait method (which performed on-demand tree traversal to derive MAP subs) was a different code path never used at this call site.

In the post-collapse path, using `fitch_subs()` is unavoidable: `collapse_edge()` calls `set_fitch_subs()` which clears `subs_ml`, so `fitch_subs()` is the only available data at that point. The `merge_shared_mutation_branches()` step runs immediately after collapse and correctly uses `fitch_subs()` for all edges. The next iteration re-runs `update_marginal()`.

## Suspected problems

- Pre-collapse decision on stale data: `find_zero_optimal_internal_edges` identifies mutation-free edges using `fitch_subs()`. An edge where Fitch assigned a substitution at an ambiguous position but marginal inference removed it would be excluded from zero-optimal candidates even though the current MAP state has no mutation on that edge. This is a conservative error (keeps edges that could be collapsed), not a correctness error (never collapses edges that carry current mutations).
- Merge comparison on mixed representations: after collapse, `merge_shared_mutation_branches()` compares `fitch_subs()` across sibling edges. For the freshly composed edge (from `chain_fitch_subs`), the Fitch representation is exact. For untouched sibling edges, `fitch_subs()` may differ from the current MAP state at ambiguous sites. This could cause missed merges (shared mutations in MAP not shared in Fitch) or incorrect merges (shared in Fitch but not in MAP). The practical impact is bounded: the next iteration re-runs marginal inference and re-identifies zero-optimal edges with the updated Fitch representation.

## Investigation tasks

- Determine whether using `subs_ml` (when available) for `find_zero_optimal_internal_edges` produces different zero-optimal edge sets on real datasets (compare edge sets from `fitch_subs()` vs `subs_ml` after `update_marginal()`).
- Determine whether using `subs_ml` for `merge_shared_mutation_branches` on non-collapsed sibling edges produces different merge decisions. This requires computing shared mutations from `subs_ml` for edges where it is available and `fitch_subs()` for freshly composed edges.
- Measure the effect on convergence speed and final topology. If Fitch-based decisions cause unnecessary iterations (conservative edge retention), switching to ML subs where available could reduce iteration count.
- Consider a hybrid approach: use `subs_ml` for `find_zero_optimal_internal_edges` (pre-collapse, all edges have marginal data), fall back to `fitch_subs()` for `merge_shared_mutation_branches` (post-collapse, composed edges lack marginal data).

## Locations

- Zero-optimal edge identification: [`packages/treetime/src/commands/optimize/run.rs#L563-L598`](../../packages/treetime/src/commands/optimize/run.rs#L563-L598) (`find_zero_optimal_internal_edges`)
- Prune and merge dispatch: [`packages/treetime/src/commands/optimize/run.rs#L599-L642`](../../packages/treetime/src/commands/optimize/run.rs#L599-L642) (`prune_and_merge_in_loop`)
- Edge collapse: [`packages/treetime/src/partition/algo/topology_cleanup/collapse.rs`](../../packages/treetime/src/partition/algo/topology_cleanup/collapse.rs)
- Shared mutation merge: [`packages/treetime/src/partition/algo/topology_cleanup/merge_shared_mutations.rs`](../../packages/treetime/src/partition/algo/topology_cleanup/merge_shared_mutations.rs) (`merge_shared_mutation_branches` and helpers)
- Fitch/ML accessors: [`packages/treetime/src/partition/sparse.rs#L112-L163`](../../packages/treetime/src/partition/sparse.rs#L112-L163) (`SparseEdgePartition`)

## Related

- [Move merge_shared_mutation_branches to shared topology_cleanup module](N-topology-cleanup-move-merge-shared-mutations.md)
- [Sparse GTR inference mixes MAP mutations with Fitch-era compositions](M-gtr-sparse-composition-stale-after-marginal.md)
