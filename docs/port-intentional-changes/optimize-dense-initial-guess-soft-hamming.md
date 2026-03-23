# Dense initial branch length guess uses soft Hamming on per-edge messages

v1 computes the dense initial branch length estimate as a soft Hamming distance using per-edge message profiles. v0 computes a soft Hamming distance using full marginal profiles evaluated at the child node. Both use the same formula (`1 - dot(pp, pc)`) but differ in the data source.

## What v0 does

v0's `optimize_tree_marginal()` ([packages/legacy/treetime/treetime/treeanc.py#L1297-L1346](../../packages/legacy/treetime/treetime/treeanc.py#L1297-L1346)) calls `marginal_branch_profile(node)` to obtain `pp` (outgroup likelihood at child) and `pc` (subtree likelihood at child), both evaluated at the child node position. The soft Hamming distance is:

```python
hamming_distance = 1 - sum(multiplicity * sum(pp * pc, axis=1)) / sum(multiplicity)
```

This distance is the bracket midpoint for Brent's method. Both profiles represent state distributions at the same node (the child), viewed from two independent information sources. Their dot product approximates the posterior certainty at each position.

## What v1 does

v1's `initial_guess_mixed()` ([packages/treetime/src/commands/optimize/optimize_unified.rs#L299-L327](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L299-L327)) calls `PartitionMarginalDense::edge_initial_differences()` ([packages/treetime/src/representation/partition/marginal_dense.rs#L124-L145](../../packages/treetime/src/representation/partition/marginal_dense.rs#L124-L145)), which computes:

```
differences = sum(1 - dot(msg_to_parent[i], msg_to_child[i]))  over non-gap positions
branch_length = differences / effective_length
```

The two messages represent state distributions at opposite ends of the edge:

- `msg_to_parent`: subtree likelihood at the child end (what the child's subtree says about the child's state)
- `msg_to_child`: outgroup likelihood at the parent end (what the rest of the tree says about the parent's state)

## Why v1 differs

v1's partition architecture stores per-edge messages rather than per-node marginal profiles. The messages are the natural data available at the point where `initial_guess_mixed` is called. Computing v0-equivalent profiles would require propagating `msg_to_child` through the transition matrix to evaluate it at the child end, adding computation without clear benefit: the initial guess is just a starting point for Newton optimization, not the final result.

The soft Hamming formula captures the key qualitative property regardless of data source: uncertain positions contribute fractionally (`1 - 1/n` for uniform profiles over `n` states) rather than as hard 0 or 1. This gives better initial estimates for Newton optimization compared to the hard argmax approach (which was the previous v1 behavior).

## Quantitative difference

For sharp profiles (one state dominates), both v0 and v1 give equivalent results: `dot ≈ 1` (agree) or `dot ≈ 0` (disagree).

For uncertain profiles, the v1 dot product at different ends of the edge and the v0 dot product at the same node give different numerical values, but both produce fractional contributions in the same range. The exact numerical difference depends on branch length and GTR model but does not affect the optimization outcome: both are reasonable starting points for Newton's method.

## Tradeoffs

**Lost:** exact numerical parity with v0's initial guess at uncertain positions.

**Gained:** the soft Hamming formula correctly handles uncertain positions, which was the main deficiency of the previous hard argmax approach. The remaining numerical difference between per-edge and per-node dot products is second-order compared to the hard-vs-soft distinction.

## Relationship to edge_subs()

`edge_initial_differences()` and `edge_subs()` serve different purposes and correctly use different data sources:

- `edge_initial_differences()` computes a continuous overlap measure (`1 - dot(pp, pc)`) for initial branch length estimation. It uses per-edge messages because the soft Hamming formula produces meaningful fractional contributions from partial-message uncertainty. This is a starting point for Newton optimization, not a final result.
- `edge_subs()` counts discrete branch mutations by comparing argmax states at parent and child nodes. It uses full node posteriors (`self.nodes[&key].profile.dis`) because discrete state assignment requires the complete marginal distribution, not a partial message whose argmax can disagree with the posterior.

The two functions must not be unified. They answer different questions with different mathematical requirements.

## Affected commands

- `optimize` (via `initial_guess_mixed`)
- `timetree` (uses `PartitionOptimizeOps` through shared optimization code)
