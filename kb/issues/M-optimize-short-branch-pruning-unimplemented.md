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

## Related material

- [kb/features/optimize.md](../features/optimize.md)
- [kb/proposals/optimize-short-branch-pruning.md](../proposals/optimize-short-branch-pruning.md)

## Related tickets

- [kb/tickets/optimize-add-short-branch-pruning.md](../tickets/optimize-add-short-branch-pruning.md)
