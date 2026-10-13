# Sparse GTR refinement counts transitions differently from dense

> [!IMPORTANT]
> **Decision required.** Which counts should sparse GTR refinement fit the model to?
>
> - Expected counts, with fixed sites weighted by their counts and one root state per site, as dense refinement does (`fn count_transitions_dense` in [packages/treetime/src/partition/marginal/shared/data.rs](../../packages/treetime/src/partition/marginal/shared/data.rs)). Sparse and dense then reach the same GTR, so dense stays an independent check of sparse. The count needs the parent and child states paired across each branch, including positions whose fixed state differs between the two sides, as `fn get_coefficients` of the sparse branch objective already does ([packages/treetime/src/partition/optimize/sparse.rs](../../packages/treetime/src/partition/optimize/sparse.rs))
> - Most likely substitutions on each branch and the composition of the most likely sequences, with the counter of the initial Fitch inference, as a commit on the `validation` branch does ([9351d4e3](https://github.com/neherlab/treetime/commit/9351d4e3c2deffe136f02a8270e790898711bfd5)). The model then fits the mutations the outputs list. The counts do not depend on the GTR, so refinement stops changing after the first iteration, and repeated substitutions at one site on one branch are not counted
>
> Evidence: the measurements below. Dense is the v0-equivalent reference: v0 counts per site (`get_branch_mutation_matrix` and the root counts of `_ml_anc_marginal` in [packages/legacy/treetime/treetime/treeanc.py](../../packages/legacy/treetime/treetime/treeanc.py)), and v1 dense starts from Fitch counts and refines with expected counts.

## Problem

GTR refinement (`--gtr-iterations` of `treetime ancestral`, and `fn refine_gtr_model_and_rate` in [packages/treetime/src/gtr/refinement.rs](../../packages/treetime/src/gtr/refinement.rs)) fits the model to expected transition counts. On the same gap-free input, the sparse and dense representations produce different counts, so the refined GTR and rate differ. Dense and sparse marginal log-likelihoods agree to about 1e-14 at the same GTR, so the difference comes from the counting step, `fn count_transitions_sparse` in [packages/treetime/src/partition/marginal/sparse/count.rs](../../packages/treetime/src/partition/marginal/sparse/count.rs).

Measured on 4 leaves, 16 gap-free sites, JC69 start, 2 iterations: dense gives mu 0.63528, sparse gives mu 0.62844.

- **Root state**: dense sums one argmax state per site (`[3, 6, 3, 4]`, sum 16). Sparse sets one one-hot vector from the argmax of the aggregated root profile (`[0, 1, 0, 0]`), so the root composition enters the equilibrium frequencies with weight 1 instead of the sequence length
- **Substitution counts**: the dense off-diagonal sum of `nij` is about 7.6, close to the 7 parsimony changes of the input; the sparse sum is about 0.57
- **Time per state**: `Ti` values are close but differ

`test_sparse_transition_counting_root_state_sums_to_length` asserts only that the root state is positive, so it does not detect the root-state difference its name describes.

## Evidence

Ignored red test: `test_refinement_sparse_gtr_and_rate_match_dense` in [packages/treetime/src/gtr/__tests__/test_refinement.rs](../../packages/treetime/src/gtr/__tests__/test_refinement.rs).

## Validation on simulated data

The validation suite (`validation/`, `ancestral` suite, variant `marginal-gtr3`: marginal reconstruction with `--gtr-iterations 3` on simulated alignments with a known GTR) shows the effect at scale: the inferred `mu` collapses to about $1/(L \mu_\text{true})$, and the equilibrium frequencies drift toward uniform. The commit on the `validation` branch names the causes in the current count: it pairs the two messages of a branch by character, so substitutions at sites fixed on both sides of the branch are never counted, and it reduces the root composition to one count.
