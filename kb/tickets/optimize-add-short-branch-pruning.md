# Add short branch pruning after optimization

After the branch-length optimization loop converges, prune internal branches too short to be resolved by the available sequence data. Matches v0's marginal-mode behavior.

## Scope

Insert between `run_optimize_loop` and the final `update_marginal` pass in `optimize/pipeline.rs`:

1. Compute `one_mutation = 1.0 / sequence_length`
2. For each internal edge: if `bl < 0.1 * one_mutation` AND P(zero) > 0.1, mark for pruning
3. Collapse marked edges (merge child into parent, creating polytomies)
4. Reconcile partition topology via `reconcile_topology`
5. Final `update_marginal` on the pruned tree

P(zero) reuses existing `is_zero_branch_optimal` from `optimize/zero_boundary.rs`.

The existing `find_zero_optimal_internal_edges` (called per iteration inside the loop) handles edges driven to exactly zero. Post-loop pruning catches edges that converged to small-but-nonzero values below the resolution of the alignment.

## Locations

- Pipeline insertion: `packages/treetime/src/optimize/pipeline.rs`
- Zero-branch logic: `packages/treetime/src/optimize/zero_boundary.rs`
- Topology collapse: `treetime-graph/src/reroot.rs` (`remove_node_if_trivial`) or existing collapse in `optimize/run_loop.rs`
- Partition reconciliation: `partition/traits.rs` (`PartitionRerootOps::reconcile_topology`)

## Tests

- flu/h3n2/20: count pruned branches, compare with v0 marginal-mode output
- Synthetic tree with short branches: verify threshold `0.1 * one_mutation`
- All branches above threshold: no pruning, identical output
- Root children below threshold: root edge must not be pruned

## Related issues

- Proposal: [kb/proposals/optimize-short-branch-pruning.md](../proposals/optimize-short-branch-pruning.md)
- Parent: [kb/proposals/optimize-pipeline-timetree-parity.md](../proposals/optimize-pipeline-timetree-parity.md)
- Feature inventory: [kb/features/optimize.md](../features/optimize.md) -- `[ ] Short branch pruning after optimization`
