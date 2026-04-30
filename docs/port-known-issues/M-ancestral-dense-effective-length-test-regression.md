# Dense/sparse optimization metric divergence with IUPAC ambiguity codes

## Failing test

`test_optimize_contribution_dense_sparse_ambiguous_r_value_and_gradient_consistency` at `packages/treetime/src/commands/optimize/__tests__/test_initial_guess_formula.rs:106`

Asserts `create_edge_contribution().evaluate(0.1)` produces identical `log_lh` and `derivative` for dense and sparse partitions at `1e-12` tolerance. Actual values: dense `-0.836`, sparse `-0.800`.

The sibling test `test_initial_guess_dense_sparse_ambiguous_r_reference_state_consistency` at line 86 (which asserts branch lengths after `initial_guess_mixed`) also fails.

## Test data

```
Tree: ((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;

A: RCGTACGT    (R = IUPAC ambiguity: A or G)
B: GCGTACGT
C: GCGTACGT
D: GCGTACGT
```

Position 0 is the only variable site. `R` in leaf A is the only ambiguity code. Both dense and sparse use the same `find_letter_ranges(seq, alphabet.unknown())` to detect unknowns, so `R` is not classified as unknown or gap in either path (it is neither `N` nor `-`).

## How each path handles R at position 0

**Sparse**: Fitch backward maps `R` to state set `{A, G}`. Intersection with B's `{G}` gives `{G}` (non-empty), so position 0 is a fixed site resolved to `G`. The sparse coefficient extraction operates on variable sites only, treating position 0 as fixed.

**Dense**: The profile map assigns `R` a mixed profile (e.g. `{A: 0.5, C: 0, G: 0.5, T: 0}`). The dense coefficient extraction operates on full N-by-K posterior matrices including position 0's mixed profile. The marginal posterior at this position incorporates evidence from the subtree (B/C/D all have G), but the mixed prior from `R` still influences the log-likelihood computation.

## Why the values diverge

The divergence is not from `edge_effective_length()` or `edge_subs()` (which are the same here since `R` is not in `non_char` for either path). It is from `create_edge_contribution()`, which computes optimization coefficients differently:

- Sparse (`packages/treetime/src/representation/partition/marginal_sparse.rs`): extracts coefficients from variable sites only. Position 0 is fixed (resolved to G by Fitch), so its contribution is a simple fixed-site term.
- Dense (`packages/treetime/src/representation/partition/marginal_dense.rs`): extracts coefficients from the full posterior matrix. Position 0 has a mixed profile from the `R` ambiguity, contributing a different log-likelihood term than a pure `G`.

The 4.5% divergence (`-0.836` vs `-0.800`) reflects the different treatment of ambiguity-code positions in the optimization objective: sparse resolves ambiguity at the Fitch stage (parsimony), dense carries it through as posterior uncertainty (marginal inference).

## Fix options

1. **Widen tolerance**: if the divergence is an acceptable consequence of dense vs sparse representation, replace `1e-12` with a tolerance that accommodates the difference. This acknowledges that dense and sparse optimization surfaces differ when ambiguity codes are present.

2. **Classify ambiguity codes as non_char in dense**: add `find_letter_ranges` calls for each ambiguity code character and merge into `non_char`. This would make dense skip ambiguous positions in `edge_subs()` and `edge_effective_length()`, matching sparse Fitch's treatment of resolved-away ambiguity sites. But it would not fix the coefficient extraction divergence.

3. **Unify coefficient extraction for ambiguous sites**: make dense treat Fitch-resolved-fixed sites the same way sparse does during coefficient extraction. This requires the dense path to distinguish "fixed by Fitch parsimony" from "fixed by sequence identity", which is not currently tracked.
