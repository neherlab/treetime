# Sparse GTR refinement counts transitions differently from dense

> [!IMPORTANT]
> **Decision required.** Which counts should sparse GTR refinement fit the model to? Both options keep Fitch counts for the initial GTR, as dense does ([kb/decisions/dense-fitch-gtr-inference.md](../decisions/dense-fitch-gtr-inference.md)):
>
> - Expected counts, with fixed sites weighted by their counts and one root state per site, as dense refinement does (`fn count_transitions_dense` [packages/treetime/src/partition/marginal/shared/data.rs#L15-L56](../../packages/treetime/src/partition/marginal/shared/data.rs#L15-L56)). Sparse and dense then reach the same GTR, so dense stays an independent check of sparse. The count needs the parent and child states paired across each branch, including positions whose fixed state differs between the two sides, as `fn get_coefficients` of the sparse branch objective already does ([packages/treetime/src/partition/optimize/sparse.rs#L11](../../packages/treetime/src/partition/optimize/sparse.rs#L11))
> - Most likely substitutions on each branch and the composition of the most likely sequences, with the counter of the initial Fitch inference. The counts do not depend on the GTR, so refinement stops changing after the first iteration, and repeated substitutions at one site on one branch are not counted. An implementation exists on the `validation` branch: [fix(treetime): count sparse GTR refinement transitions from branch mutations](https://github.com/neherlab/treetime/commit/9351d4e3c2deffe136f02a8270e790898711bfd5)
>
> Evidence: the measurements below. Dense is the v0-equivalent reference: v0 counts per site (`def get_branch_mutation_matrix` [packages/legacy/treetime/treetime/treeanc.py#L1085](../../packages/legacy/treetime/treetime/treeanc.py#L1085) and the root counts of `def _ml_anc_marginal` [packages/legacy/treetime/treetime/treeanc.py#L762](../../packages/legacy/treetime/treetime/treeanc.py#L762)), and v1 dense starts from Fitch counts and refines with expected counts.

## Problem

GTR refinement (`--gtr-iterations` of `treetime ancestral`, and `fn refine_gtr_model_and_rate` in [packages/treetime/src/gtr/refinement.rs](../../packages/treetime/src/gtr/refinement.rs)) fits the model to expected transition counts. On the same gap-free input, the sparse and dense representations produce different counts, so the refined GTR and rate differ. Dense and sparse marginal log-likelihoods agree to about 1e-14 at the same GTR, so the difference comes from the counting step, `fn count_transitions_sparse` in [packages/treetime/src/partition/marginal/sparse/count.rs](../../packages/treetime/src/partition/marginal/sparse/count.rs).

Measured on 4 leaves, 16 gap-free sites, JC69 start, 2 iterations: dense gives mu 0.63528, sparse gives mu 0.62844.

- **Root state**: dense sums one argmax state per site (`[3, 6, 3, 4]`, sum 16). Sparse sets one one-hot vector from the argmax of the aggregated root profile (`[0, 1, 0, 0]`), so the root composition enters the equilibrium frequencies with weight 1 instead of the sequence length
- **Substitution counts**: the dense off-diagonal sum of `nij` is about 7.6, close to the 7 parsimony changes of the input; the sparse sum is about 0.57
- **Time per state**: `Ti` values are close but differ

`test_sparse_transition_counting_root_state_sums_to_length` asserts only that the root state is positive, so it does not detect the root-state difference its name describes.

## Evidence

On simulated alignments with a known GTR and `--gtr-iterations 3` (variant `marginal-gtr3` of the ancestral suite in [validation/README.md](../../validation/README.md)), with alignment length $L$ and simulated substitution rate $\mu_{\text{true}}$, the inferred sparse rate collapses to about $1/(L \mu_{\text{true}})$ and the equilibrium frequencies drift toward uniform. Two defects of `fn count_transitions_sparse` [packages/treetime/src/partition/marginal/sparse/count.rs#L16](../../packages/treetime/src/partition/marginal/sparse/count.rs#L16) explain it: it pairs the two messages of a branch by character, so it never counts a substitution at a site that is fixed on both sides, and it reduces the root composition to one count.

Ignored red test: `test_refinement_sparse_gtr_and_rate_match_dense` in [packages/treetime/src/gtr/__tests__/test_refinement.rs](../../packages/treetime/src/gtr/__tests__/test_refinement.rs).

## Ambiguous leaf positions count twice

`fn accumulate_sparse_transitions()` [packages/treetime/src/partition/marginal/sparse/count.rs#L58-L85](../../packages/treetime/src/partition/marginal/sparse/count.rs#L58-L85) sums the expected transitions of an edge over the explicit positions of `msg_to_parent` with weight 1, and over its fixed states with the weight of each fixed `count`.

On a leaf edge, `msg_to_parent` is the leaf message (`fn leaf_message()` [packages/treetime/src/partition/marginal/sparse/backward.rs#L96-L119](../../packages/treetime/src/partition/marginal/sparse/backward.rs#L96-L119)). Its explicit positions are the ambiguous positions of the leaf. Its fixed profiles cover every determined state, ambiguous states included, each with the composition count of that state in the leaf sequence, because `fn observe_seq()` puts the ambiguous positions into `variable` and also counts every character into the composition [packages/treetime/src/partition/storage/sparse.rs#L285-L316](../../packages/treetime/src/partition/storage/sparse.rs#L285-L316). An `R` position of a leaf therefore enters twice: once as an explicit position with the profile of `R`, and once in the fixed group of `R`, whose count includes it. Internal messages do not have this overlap, because their fixed counts exclude the positions that stay explicit.

Both options of the decision above rewrite this counter, so the fix of the double count is part of the same change. No test covers an ambiguous leaf position yet.

## Validation

- `test_refinement_sparse_gtr_and_rate_match_dense` passes and loses its `#[ignore]`
- A tree with one ambiguous leaf position gives transition counts `nij` and times per state `Ti` equal to a hand computation of the expected counts
- `test_sparse_transition_counting_root_state_sums_to_length` asserts the sum it names
