# Dense initial branch length guess uses hard substitution count

v1 computes the dense initial branch length estimate as a hard substitution count divided by effective alignment length. v0 computes a soft Hamming distance using full marginal profiles evaluated at the child node. The two approaches produce the same result for sharp profiles but differ for uncertain positions.

## What v0 does

v0's `optimize_tree_marginal()` ([packages/legacy/treetime/treetime/treeanc.py#L1297-L1346](../../packages/legacy/treetime/treetime/treeanc.py#L1297-L1346)) calls `marginal_branch_profile(node)` to obtain `pp` (outgroup likelihood at child) and `pc` (subtree likelihood at child), both evaluated at the child node position. The soft Hamming distance is:

```python
hamming_distance = 1 - sum(multiplicity * sum(pp * pc, axis=1)) / sum(multiplicity)
```

This distance is the bracket midpoint for Brent's method. Both profiles represent state distributions at the same node (the child), viewed from two independent information sources. Their dot product approximates the posterior certainty at each position.

## What v1 does

v1's `initial_guess_mixed()` ([packages/treetime/src/commands/optimize/optimize_unified.rs#L729-L806](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L729-L806)) computes:

```
sub_count    = sum of edge_subs().len() across partitions
eff_length   = sum of edge_effective_length() across partitions
branch_length = sub_count / eff_length
```

`edge_subs()` on dense partitions ([packages/treetime/src/representation/partition/marginal_dense.rs#L81-L112](../../packages/treetime/src/representation/partition/marginal_dense.rs#L81-L112)) takes the hard `argmax` of parent and child profiles and counts positions where the MAP states differ (excluding gaps, ambiguous, and non-canonical states). `edge_effective_length()` counts non-gap positions.

This is a hard Hamming distance: each position contributes 0 (same MAP state) or 1 (different MAP state), with no fractional contributions from uncertain profiles.

## Why v1 differs

v1's partition architecture stores per-edge messages rather than per-node marginal profiles. Computing v0-equivalent dot products would require propagating messages through the transition matrix to evaluate both at the same node position.

The initial guess serves only as a starting point for Newton/Brent optimization. Hard counting is adequate because:

- For sharp profiles (one state dominates), hard and soft approaches agree exactly
- For uncertain profiles, the hard count slightly overestimates distance (uncertain positions that happen to have different argmax states contribute 1.0 instead of a fractional value), but this is a conservative starting point that Newton optimization corrects

## Tradeoffs

**Lost:** fractional distance contributions from uncertain positions. v0's soft Hamming gives a smoother, more accurate initial estimate when profiles are diffuse.

**Gained:** simpler implementation that reuses `edge_subs()` (shared between sparse and dense partitions). No sensitivity to profile shape at uncertain positions.

## Affected commands

- `optimize` (via `initial_guess_mixed`)
- `timetree` (uses `PartitionOptimizeOps` through shared optimization code)
