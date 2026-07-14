# Add short branch pruning after optimization

After the branch-length optimization loop converges, prune internal branches too short to be resolved by the available sequence data. Matches v0's marginal-mode behavior.

## Scope

Insert between `run_optimize_loop` and the final `update_marginal` pass in `optimize/pipeline.rs`:

1. Let $L$ be the sequence length and compute the one-mutation resolution $m=1/L$.
2. For each eligible internal edge, mark it when its branch length $b$ satisfies $b < 0.1m$ and the zero-time sequence transition likelihood of its reconstructed parent and child sequences satisfies $P(0) > 0.1$.
3. Collapse marked edges (merge child into parent, creating polytomies)
4. Reconcile partition topology via `reconcile_topology`
5. Final `update_marginal` on the pruned tree

Compute $P(0)$ with the same sequence states and pattern multiplicities as `def TreeAnc.prune_short_branches()`. [`packages/legacy/treetime/treetime/treeanc.py#L1475-L1495`](../../packages/legacy/treetime/treetime/treeanc.py#L1475-L1495) `def TreeAnc.optimize_tree()` supplies the marginal-mode orchestration. [`packages/legacy/treetime/treetime/treeanc.py#L1424-L1444`](../../packages/legacy/treetime/treetime/treeanc.py#L1424-L1444) Do not substitute `fn is_zero_branch_optimal()`: it tests a likelihood derivative and has a different model domain.

The existing `find_zero_optimal_internal_edges` (called per iteration inside the loop) handles edges driven to exactly zero. Post-loop pruning catches edges that converged to small-but-nonzero values below the resolution of the alignment.

## Locations

- Pipeline insertion: `packages/treetime/src/optimize/pipeline.rs`
- Zero-branch logic: `packages/treetime/src/optimize/zero_boundary.rs`
- Topology collapse: `fn collapse_edge()` [`packages/treetime/src/optimize/topology/collapse.rs#L34-L82`](../../packages/treetime/src/optimize/topology/collapse.rs#L34-L82)
- Partition reconciliation: `partition/traits.rs` (`PartitionRerootOps::reconcile_topology`)

## Tests

- flu/h3n2/20: count pruned branches, compare with v0 marginal-mode output
- Synthetic tree with short branches: verify both the $0.1m$ length threshold and the exact $P(0)$ threshold.
- Cover JC69 and at least one model for which `is_zero_branch_optimal()` declines to decide.
- Include a fixture where the probability threshold and derivative-sign predicate disagree.
- All branches above threshold: no pruning, identical output
- Root has no incoming edge and is never a pruning candidate.
- Terminal children are excluded.
- Eligible internal children of the root use the same predicate as other internal nodes.

## Related issues

- Source: [kb/issues/M-optimize-short-branch-pruning-unimplemented.md](../issues/M-optimize-short-branch-pruning-unimplemented.md) -- delete after full resolution
- Proposal: [kb/proposals/optimize-short-branch-pruning.md](../proposals/optimize-short-branch-pruning.md)
- Parent: [kb/proposals/optimize-pipeline-timetree-parity.md](../proposals/optimize-pipeline-timetree-parity.md)
- Feature inventory: [kb/features/optimize.md](../features/optimize.md) -- `[ ] Short branch pruning after optimization`
