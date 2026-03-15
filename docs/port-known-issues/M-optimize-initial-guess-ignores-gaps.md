# Initial branch length guess ignores gaps and unknowns

The optimize command's initial branch length estimate counts raw substitution differences without excluding gap or ambiguous positions, inflating branch lengths for gappy alignments.

## Problem

`initial_guess_mixed()` ([packages/treetime/src/commands/optimize/optimize_unified.rs#L257-L280](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L257-L280)) computes the initial branch length as:

```
branch_length = differences / total_length
```

Where `differences` counts substitutions from `edge_subs()` (Fitch parsimony output for sparse, sequence comparison for dense) and `total_length` sums partition sequence lengths. Neither value filters gap or ambiguous positions.

A gap-to-nucleotide pair registers as a substitution, but it represents missing data, not evolutionary distance. An alignment with 10% gaps can inflate the initial branch length estimate by up to 10%, pushing Newton's method into a region where convergence is slower or the Hessian flips sign, triggering the grid search fallback.

## v0 handling

v0 computes Hamming distance in `optimal_t_compressed()` via `state_pair()` ([packages/legacy/treetime/treetime/gtr.py#L631-L712](../../packages/legacy/treetime/treetime/gtr.py#L631-L712)), which accepts an `ignore_gaps` parameter. When `ignore_gaps=True` (the default for branch length optimization), gap positions are excluded from both the numerator (substitution count) and denominator (total pair count).

v0 also applies `MIN_BRANCH_LENGTH` floor (`1e-3 * one_mutation`) in `_branch_length_to_gtr()` ([packages/legacy/treetime/treetime/treeanc.py#L752-L760](../../packages/legacy/treetime/treetime/treeanc.py#L752-L760)) before any GTR calculation. This prevents numerically degenerate branch lengths (see also `N-core-branch-length-clamping.md`).

## Proposed solutions

### S1: Filter gaps in initial guess (minimal)

In `initial_guess_mixed()`, exclude gap-to-anything substitutions from the `differences` count. For sparse partitions, filter `edge_subs()` results where either the reference or query state maps to a gap character. For dense partitions, reduce `total_length` by the number of gap positions per edge.

### S2: Use effective alignment length (full fix)

Compute per-edge effective alignment length: the number of positions where both parent and child have non-gap, non-ambiguous states. Use this as the denominator instead of the full partition length. This matches v0's semantics where `multiplicity` sums exclude masked positions.

## Interaction with other issues

- [Optimize sparse marginal crashes with NaN](H-optimize-sparse-marginal-nan-crash.md): inflated initial branch lengths can push the marginal reconstruction into underflow territory on early iterations, contributing to the NaN cascade
- [Gap character not handled in alphabet](H-timetree-gap-alphabet.md): related gap handling inconsistencies across commands
