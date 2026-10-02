# Sparse GTR refinement counts transitions differently from dense

## Problem

GTR refinement (`--gtr-iterations` of `treetime ancestral`, and `fn refine_gtr_model_and_rate` in [packages/treetime/src/gtr/refinement.rs](../../packages/treetime/src/gtr/refinement.rs)) fits the model to expected transition counts. On the same gap-free input, the sparse and dense representations produce different counts, so the refined GTR and rate differ. Dense and sparse marginal log-likelihoods agree to about 1e-14 at the same GTR, so the difference comes from the counting step, `fn count_transitions_sparse` in [packages/treetime/src/partition/marginal/sparse/count.rs](../../packages/treetime/src/partition/marginal/sparse/count.rs).

Measured on 4 leaves, 16 gap-free sites, JC69 start, 2 iterations: dense gives mu 0.63528, sparse gives mu 0.62844.

- **Root state**: dense sums one argmax state per site (`[3, 6, 3, 4]`, sum 16). Sparse sets one one-hot vector from the argmax of the aggregated root profile (`[0, 1, 0, 0]`), so the root composition enters the equilibrium frequencies with weight 1 instead of the sequence length
- **Substitution counts**: the dense off-diagonal sum of `nij` is about 7.6, close to the 7 parsimony changes of the input; the sparse sum is about 0.57
- **Time per state**: `Ti` values are close but differ

`test_sparse_transition_counting_root_state_sums_to_length` asserts only that the root state is positive, so it does not detect the root-state difference its name describes.

## Evidence

Ignored red test: `test_refinement_sparse_gtr_and_rate_match_dense` in [packages/treetime/src/gtr/__tests__/test_refinement.rs](../../packages/treetime/src/gtr/__tests__/test_refinement.rs).

## Open question

Dense is the v0-equivalent reference (v0 counts per site, `packages/legacy/treetime/treetime/treeanc.py`, `get_branch_mutation_matrix` and `_ml_anc_marginal` root counts). Decide whether sparse must reproduce the dense counts exactly, including fixed sites weighted by their counts and one root state per site, and then fix `count_transitions_sparse`.
