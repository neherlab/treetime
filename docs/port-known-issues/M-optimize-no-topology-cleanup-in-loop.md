# Optimize loop does not prune zero-length branches or resolve polytomies

v0's `optimize_tree()` calls `prune_short_branches()` inside each joint-mode iteration, collapsing internal edges whose optimal length is zero. v1's `run_optimize()` detects zero-optimal branches via `is_zero_branch_optimal()` and sets `branch_length = 0.0`, but the edge remains in the tree. No topology changes occur during optimization.

## v0 behavior

`optimize_tree()` ([`packages/legacy/treetime/treetime/treeanc.py#L1384-L1473`](../../packages/legacy/treetime/treetime/treeanc.py#L1384-L1473)) with `prune_short=True` (default):

- Joint mode: `prune_short_branches()` called inside each iteration (line 1458-1459), before ancestral reconstruction
- Marginal mode: `prune_short_branches()` called once after `optimize_tree_marginal()` completes (line 1436)

`prune_short_branches()` ([`packages/legacy/treetime/treetime/treeanc.py#L1475-L1496`](../../packages/legacy/treetime/treetime/treeanc.py#L1475-L1496)) uses a compound criterion:

1. `branch_length < 0.1 * one_mutation` (less than 10% of one expected mutation)
2. `gtr.prob_t(parent_seq, child_seq, 0.0) > 0.1` (likelihood-based check that sequences are compatible at zero distance)

Both conditions must hold. Condition 2 prevents collapsing branches that are short but carry genuine signal.

## v1 behavior

`run_optimize()` ([`packages/treetime/src/commands/optimize/run.rs#L127-L148`](../../packages/treetime/src/commands/optimize/run.rs#L127-L148)):

- `is_zero_branch_optimal()` ([`packages/treetime/src/commands/optimize/optimize_unified.rs#L192-L214`](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L192-L214)) detects zero-optimal branches using the derivative-at-zero sign
- When detected, sets `branch_length = 0.0` but does NOT collapse the edge
- Zero-length edges accumulate across iterations
- No polytomy resolution after pruning (design doc [`docs/algorithms/optimize.md#L19-L20`](../algorithms/optimize.md#L19-L20) describes both as "added value")

## Impact

- Wasted computation on degenerate zero-length edges in subsequent iterations
- No polytomy formation, so shared-mutation merging (the complementary operation) cannot run
- Optimization converges differently from v0: v0's per-iteration pruning simplifies the tree progressively

## Blocker

`M-prune-collapse-uses-union-not-composition.md`: edge collapsing uses substitution union instead of composition. Must be fixed before integrating edge collapsing into the optimize loop.

## Design doc reference

[`docs/algorithms/optimize.md#L19-L20`](../algorithms/optimize.md#L19-L20):

- "prune zero length branches"
- "merge branches in polytomies that share mutations"

See [`docs/reports/iterative-tree-refinement/_index.md`](../reports/iterative-tree-refinement/_index.md) for full analysis (chapters 6, 9, and 10).

Supersedes previous `N-optimize-polytomy-branch-merging` (Negligible) -- reclassified as Medium because the missing topology cleanup affects optimization convergence and downstream operations.
