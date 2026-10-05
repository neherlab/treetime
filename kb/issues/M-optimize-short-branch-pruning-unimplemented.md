# Optimize omits v0 post-convergence short-branch pruning

The optimize pipeline does not prune internal branches whose small nonzero lengths cannot be resolved by the available sequence data. TreeTime v0 marginal mode performs this cleanup after branch-length optimization; the v1 feature inventory records it as unimplemented.

## Reference behavior

`def TreeAnc.optimize_tree()` invokes `def TreeAnc.prune_short_branches()` after marginal branch optimization. [`packages/legacy/treetime/treetime/treeanc.py#L1424-L1444`](../../packages/legacy/treetime/treetime/treeanc.py#L1424-L1444) The pruning function removes an eligible internal edge when both conditions hold. [`packages/legacy/treetime/treetime/treeanc.py#L1475-L1495`](../../packages/legacy/treetime/treetime/treeanc.py#L1475-L1495)

- let $L$ be the sequence length and $m=1/L$ the one-mutation resolution; the branch length $b$ satisfies $b < 0.1m$; and
- the zero-time sequence transition likelihood for the reconstructed parent and child sequences satisfies $P(0) > 0.1$.

V1 already identifies some zero-optimal edges during optimization through `is_zero_branch_optimal()`. That derivative-sign test is not equivalent to v0's probability threshold and deliberately declines to decide for several substitution models. The missing post-convergence cleanup therefore requires the exact v0 probability calculation rather than reuse of the existing predicate.

## Impact

Unresolved internal edges remain in optimized trees instead of being collapsed into polytomies. This diverges from v0 marginal-mode topology cleanup and can retain branch lengths below the alignment's resolution.

## Expected behavior

After `run_optimize_loop()` and before the final marginal update:

1. identify non-root internal edges satisfying both v0 pruning conditions;
2. collapse the marked edges;
3. reconcile partition topology; and
4. recompute the final marginal state on the pruned tree.

Compute $P(0)$ with the same sequence states and pattern multiplicities as `def TreeAnc.prune_short_branches()`. Do not substitute `fn is_zero_branch_optimal()`: it tests a likelihood derivative and has a different model domain.

The existing `fn find_zero_optimal_internal_edges()` [packages/treetime/src/optimize/run_loop.rs#L236](../../packages/treetime/src/optimize/run_loop.rs#L236), called in every iteration of the loop, handles edges driven to exactly zero. Post-loop pruning catches edges that converged to small but nonzero values below the resolution of the alignment.

## Locations

- Pipeline insertion: `run_optimize_loop()` call and final `marginal_update()` [packages/treetime/src/optimize/pipeline.rs#L136-L155](../../packages/treetime/src/optimize/pipeline.rs#L136-L155)
- Zero-branch predicate: `fn is_zero_branch_optimal()` [packages/treetime/src/optimize/zero_boundary.rs#L16](../../packages/treetime/src/optimize/zero_boundary.rs#L16)
- Topology collapse: `fn collapse_edge()` [packages/treetime/src/optimize/topology/collapse.rs#L7-L43](../../packages/treetime/src/optimize/topology/collapse.rs#L7-L43)
- Partition reconciliation: `fn MarginalReconstruction::reconcile_topology()` [packages/treetime/src/partition/marginal/reconstruction.rs#L286](../../packages/treetime/src/partition/marginal/reconstruction.rs#L286)

## Validation

- flu/h3n2/20: count pruned branches and compare with v0 marginal-mode output
- Synthetic tree with short branches: verify both the $0.1m$ length threshold and the exact $P(0)$ threshold
- Cover JC69 and at least one model for which `is_zero_branch_optimal()` declines to decide
- Include a fixture where the probability threshold and the derivative-sign predicate disagree
- All branches above threshold: no pruning, identical output
- The root has no incoming edge and is never a pruning candidate
- Terminal children are excluded
- Eligible internal children of the root use the same predicate as other internal nodes

## Related material

- [kb/features/optimize.md](../features/optimize.md) -- `[ ] Short branch pruning after optimization`
- [kb/proposals/optimize-short-branch-pruning.md](../proposals/optimize-short-branch-pruning.md)
- [kb/proposals/optimize-pipeline-timetree-parity.md](../proposals/optimize-pipeline-timetree-parity.md)
