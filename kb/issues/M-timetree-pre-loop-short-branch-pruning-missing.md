# Timetree omits v0 pre-loop short-branch pruning

The v1 timetree pipeline never collapses short internal branches. TreeTime v0 collapses them once before the timetree iterations, so that arbitrary bifurcations from tree builders become polytomies before temporal inference and polytomy resolution.

## Reference behavior

Tree builders such as IQ-TREE, FastTree, and RAxML resolve unsupported splits into bifurcations with zero or near-zero internal branch lengths. The order of these bifurcations is arbitrary and is usually inconsistent with the sampling dates. v0 removes such branches before the timetree loop:

- Input branch lengths: `def TreeTime._run()` calls `def TreeAnc.prune_short_branches()` after ancestral reconstruction [`packages/legacy/treetime/treetime/treetime.py#L235-L241`](../../packages/legacy/treetime/treetime/treetime.py#L235-L241)
- Optimized branch lengths: `def TreeTime._run()` calls `optimize_tree(max_iter=1)` [`packages/legacy/treetime/treetime/treetime.py#L242-L243`](../../packages/legacy/treetime/treetime/treetime.py#L242-L243), which prunes after optimization because `prune_short` defaults to `True` [`packages/legacy/treetime/treetime/treetime.py#L215`](../../packages/legacy/treetime/treetime/treetime.py#L215), [`packages/legacy/treetime/treetime/treeanc.py#L1432-L1444`](../../packages/legacy/treetime/treetime/treeanc.py#L1432-L1444)
- After rerooting: the second pre-loop `optimize_tree(max_iter=1)` call prunes again [`packages/legacy/treetime/treetime/treetime.py#L262-L266`](../../packages/legacy/treetime/treetime/treetime.py#L262-L266)
- Inside the loop: after polytomy resolution v0 sets `prune_short = False` [`packages/legacy/treetime/treetime/treetime.py#L324-L329`](../../packages/legacy/treetime/treetime/treetime.py#L324-L329), so later calls keep the resolved topology

`def TreeAnc.prune_short_branches()` [`packages/legacy/treetime/treetime/treeanc.py#L1475-L1495`](../../packages/legacy/treetime/treetime/treeanc.py#L1475-L1495) removes a non-root internal edge when both conditions hold, with $L$ the sequence length and $m = 1/L$ the one-mutation resolution:

- the branch length $b$ satisfies $b < 0.1m$
- the probability of the reconstructed parent and child sequences at zero time, with pattern multiplicities, satisfies $P(0) > 0.1$

## Current behavior

The v1 pre-loop (`fn run_pre_loop()` [`packages/treetime/src/timetree/pre_loop.rs#L36`](../../packages/treetime/src/timetree/pre_loop.rs#L36)) optimizes branch lengths (`fn ml_optimize()`, `fn optimize_branch_lengths()`), reroots, and filters clock outliers. None of these steps collapses an edge, and no other code in `packages/treetime/src/timetree/` collapses or prunes edges. [kb/decisions/timetree-no-zero-branch-collapse-in-loop.md](../decisions/timetree-no-zero-branch-collapse-in-loop.md) excludes collapse inside the loop only.

## Impact

Arbitrary bifurcations from the input tree stay in the topology. Temporal polytomy resolution (`--resolve-polytomies`) acts only on existing polytomies, so it cannot reorder these splits by the sampling dates. On trees with many near-identical sequences, v1 timetrees can therefore keep clades that conflict with the temporal order and differ in topology from v0 output.

## Expected behavior

Before the timetree loop, after the pre-loop ancestral reconstruction or branch-length optimization:

1. identify non-root internal edges satisfying both v0 conditions
2. collapse the marked edges
3. reconcile partition topology
4. recompute the marginal state on the pruned tree

The criterion is the same as in [M-optimize-short-branch-pruning-unimplemented.md](M-optimize-short-branch-pruning-unimplemented.md); both commands should share one implementation.

## Locations

- Pre-loop: `fn run_pre_loop()` [`packages/treetime/src/timetree/pre_loop.rs#L36`](../../packages/treetime/src/timetree/pre_loop.rs#L36)
- Topology collapse: `fn collapse_edge()` [`packages/treetime/src/optimize/topology/collapse.rs`](../../packages/treetime/src/optimize/topology/collapse.rs)

## Validation

- flu/h3n2/20 and a dataset with many identical sequences: count pruned branches and compare topology with v0 `treetime` output with default options
