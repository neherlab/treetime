# Dense effective length change breaks two tests

After adding `unknown`/`non_char` tracking to dense partitions, `edge_effective_length()` subtracts `non_char` (gaps + unknowns) instead of gaps only. This is the correct behavior (matches sparse), but two pre-existing tests have hardcoded expectations derived from the old gaps-only computation.

## Failing tests

1. `test_marginal_dense_update_is_idempotent` at `packages/treetime/src/commands/ancestral/__tests__/test_marginal_dense.rs:312`: golden log-likelihood `-57.712498930787206` was computed with gaps-only effective length. The corrected effective length shifts it to approximately `-57.83`.

2. `test_optimize_contribution_dense_sparse_ambiguous_r_value_and_gradient_consistency` at `packages/treetime/src/commands/optimize/__tests__/test_initial_guess_formula.rs:125`: dense and sparse log-likelihoods diverge at `1e-12` tolerance. The sparse path already subtracted `non_char`; dense previously subtracted gaps only. The correction brings them closer in principle but the test data contains IUPAC ambiguity codes (`R`) that interact differently with `non_char` tracking in dense vs sparse representations.

## Root cause

`edge_effective_length()` at `packages/treetime/src/representation/partition/marginal_dense.rs:116` changed from subtracting gap positions to subtracting `non_char` positions. For sequences without unknowns, `non_char == gaps` and the result is identical. For sequences with IUPAC ambiguity codes or `N` characters, `non_char` is a superset of gaps, reducing effective length. The backward-pass `non_char` computation (intersection of children's `non_char`) can also produce different ranges than the backward-pass `gaps` computation when children disagree on gap vs unknown status.

## Fix

1. Update the golden log-likelihood in `test_marginal_dense_update_is_idempotent` to match the corrected value. Verify the new value against v0 oracle.
2. Investigate the dense/sparse divergence in the ambiguous-R test. Determine whether the remaining divergence is within expected numerical tolerance or indicates a gap/unknown classification discrepancy for IUPAC codes.
