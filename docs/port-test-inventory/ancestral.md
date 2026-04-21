# Ancestral Reconstruction Tests

[Back to index](_index.md)

## Summary

| Category                 | Files  | Tests   | Type            |
| ------------------------ | ------ | ------- | --------------- |
| Fitch parsimony          | 1      | 7       | Unit            |
| Marginal ML dense        | 1      | 5       | Unit            |
| Marginal ML sparse       | 1      | 6       | Unit            |
| Dense/sparse equivalence | 2      | 2       | Unit + Property |
| Idempotency              | 2      | 4       | Unit + Property |
| Normalization            | 2      | 6       | Unit + Property |
| Root invariance          | 1      | 4       | Property + Unit |
| Python parity            | 1      | 8       | Unit            |
| Consistency              | 1      | 5       | Unit            |
| Branch length            | 2      | 12      | Unit            |
| Topology                 | 3      | 13      | Unit            |
| Stability                | 3      | 16      | Unit            |
| Analytical               | 3      | 10      | Golden-master   |
| Sparse composition       | 1      | 12      | Unit            |
| Generator validation     | 4      | 15      | Property        |
| **Total**                | **28** | **125** | Mixed           |

Support files (helpers only, no tests): `prop_marginal_support.rs`, `test_marginal_analytical_support.rs`, `test_marginal_stability_support.rs`

---

## Fitch Parsimony

**File:** [`test_fitch.rs`](../../packages/treetime/src/commands/ancestral/__tests__/test_fitch.rs)

| Test                                                       | Purpose                                                     |
| ---------------------------------------------------------- | ----------------------------------------------------------- |
| `test_ancestral_reconstruction_fitch`                      | MAP sequences at internal nodes                             |
| `test_ancestral_reconstruction_fitch_with_leaves`          | Same with leaf sequences included                           |
| `test_compress_sequences_retains_internal_exact_sequences` | Internal exact sequences remain available after compression |
| `test_fitch_internals`                                     | Substitution and indel mutations on edges                   |
| `test_fitch_complex_gaps`                                  | Overlapping deletions and variable insertions               |
| `test_fitch_polytomy`                                      | Multifurcation with 3 children                              |
| `test_fitch_backward_state`                                | Intermediate state between backward/forward passes          |

**Algorithm:** Fitch maximum parsimony (backward pass, forward pass, gap tracking)

---

## Marginal ML - Dense

**File:** [`test_marginal_dense.rs`](../../packages/treetime/src/commands/ancestral/__tests__/test_marginal_dense.rs)

| Test                                                          | Purpose                                            |
| ------------------------------------------------------------- | -------------------------------------------------- |
| `test_ancestral_reconstruction_marginal_dense`                | MAP sequences at internal nodes under JC69         |
| `test_marginal_dense_probability_normalization`               | All posteriors sum to 1.0 within 4 ULPs            |
| `test_marginal_dense_update_is_idempotent`                    | Second update_marginal returns same log-likelihood |
| `test_marginal_dense_log_lh_root_invariance_reversible_model` | Pulley principle across three rootings             |
| `test_total_likelihood_marginal_dense_all_triplets`           | Law of total probability over 64 triplets          |

**Algorithm:** Felsenstein's pruning (sum-product belief propagation on trees)

---

## Marginal ML - Sparse

**File:** [`test_marginal_sparse.rs`](../../packages/treetime/src/commands/ancestral/__tests__/test_marginal_sparse.rs)

| Test                                                           | Purpose                                           |
| -------------------------------------------------------------- | ------------------------------------------------- |
| `test_ancestral_reconstruction_marginal_sparse`                | MAP sequences using Fitch compression + marginal  |
| `test_marginal_sparse_probability_normalization`               | Variable and fixed distributions sum to 1.0       |
| `test_marginal_sparse_update_is_idempotent`                    | Fixed-point test for sparse representation        |
| `test_marginal_sparse_log_lh_root_invariance_reversible_model` | Root invariance with non-uniform GTR              |
| `test_marginal_sparse_posterior_values_python_parity`          | Specific posterior values match Python v0         |
| `test_total_likelihood_marginal_sparse_all_triplets`           | Law of total probability for sparse (64 triplets) |

**Algorithm:** Fitch compression + marginal ML on variable positions

---

## Dense/Sparse Equivalence

**Files:** [`test_marginal_dense_sparse_example.rs`](../../packages/treetime/src/commands/ancestral/__tests__/test_marginal_dense_sparse_example.rs), [`test_marginal_dense_sparse_prop.rs`](../../packages/treetime/src/commands/ancestral/__tests__/test_marginal_dense_sparse_prop.rs)

| Test                                                      | Purpose                                               | Type     |
| --------------------------------------------------------- | ----------------------------------------------------- | -------- |
| `test_marginal_dense_sparse_example_gap_free_consistency` | Dense and sparse log-likelihoods agree on fixed input | Unit     |
| `test_prop_marginal_dense_sparse_gap_free_consistency`    | Dense and sparse agree across 30 random inputs        | Property |

**Purpose:** Cross-validate sparse optimization against dense reference
**Known issue:** ~2.5% of random GTR configs show divergence up to 7.6e-6 relative difference

---

## Idempotency Tests

**Files:** [`test_marginal_idempotency_example.rs`](../../packages/treetime/src/commands/ancestral/__tests__/test_marginal_idempotency_example.rs), [`test_marginal_idempotency_prop.rs`](../../packages/treetime/src/commands/ancestral/__tests__/test_marginal_idempotency_prop.rs)

| Test                                       | Purpose                                                           | Type     |
| ------------------------------------------ | ----------------------------------------------------------------- | -------- |
| `test_marginal_idempotency_example_dense`  | Fixed tree: update_marginal twice yields identical log-likelihood | Unit     |
| `test_marginal_idempotency_example_sparse` | Same for sparse                                                   | Unit     |
| `test_prop_marginal_idempotency_dense`     | Idempotency across 50 random inputs                               | Property |
| `test_prop_marginal_idempotency_sparse`    | Idempotency for sparse across 50 random inputs                    | Property |

**Invariant:** Sum-product on trees converges in one pass

---

## Normalization Tests

**Files:** [`test_marginal_normalization_example.rs`](../../packages/treetime/src/commands/ancestral/__tests__/test_marginal_normalization_example.rs), [`test_marginal_normalization_prop.rs`](../../packages/treetime/src/commands/ancestral/__tests__/test_marginal_normalization_prop.rs)

| Test                                                    | Purpose                                          | Type     |
| ------------------------------------------------------- | ------------------------------------------------ | -------- |
| `test_marginal_normalization_example_dense`             | All posteriors sum to 1.0, finite, non-negative  | Unit     |
| `test_marginal_normalization_example_sparse`            | Same for sparse (variable + fixed distributions) | Unit     |
| `test_prop_marginal_normalization_dense`                | Normalization across 50 small random inputs      | Property |
| `test_prop_marginal_normalization_sparse`               | Sparse normalization across 50 small inputs      | Property |
| `test_prop_marginal_normalization_dense_log_lh_finite`  | Log-likelihood finite and <=0 for larger inputs  | Property |
| `test_prop_marginal_normalization_sparse_log_lh_finite` | Same for sparse                                  | Property |

**Invariant:** P(s|D) is a valid probability distribution at every node and position

---

## Root Invariance Tests

**File:** [`test_marginal_root_invariance_prop.rs`](../../packages/treetime/src/commands/ancestral/__tests__/test_marginal_root_invariance_prop.rs)

| Test                                                     | Purpose                                                 | Type     |
| -------------------------------------------------------- | ------------------------------------------------------- | -------- |
| `test_prop_marginal_dense_log_lh_root_invariance`        | Log-likelihood invariant under rerooting (30 cases)     | Property |
| `test_prop_marginal_sparse_log_lh_root_invariance`       | Same for sparse (tolerance 1e-1 due to Fitch ambiguity) | Property |
| `test_reroot_at_internal_node_preserves_topology`        | Rerooting helper preserves leaves and topology          | Unit     |
| `test_reroot_at_internal_node_different_rootings_differ` | Different indices produce different rootings            | Unit     |

**Invariant:** Felsenstein's pulley principle for reversible models. Helper tests are in `mod helpers::tests` within the same file.

---

## Python v0 Parity Tests

**File:** [`test_python_parity.rs`](../../packages/treetime/src/commands/ancestral/__tests__/test_python_parity.rs)

| Test                                            | Purpose                                                        |
| ----------------------------------------------- | -------------------------------------------------------------- |
| `test_root_sequence_matches_python_h3n2_na_20`  | Root sequence matches Python v0 on H3N2 dataset                |
| `test_internal_node_ab_profile_matches_python`  | Node AB position 0 posterior matches Python                    |
| `test_root_profile_matches_python`              | Root position 0 posterior matches Python                       |
| `test_internal_node_cd_profile_valid`           | Node CD posteriors normalized                                  |
| `test_all_internal_nodes_normalized`            | All internal node posteriors normalized                        |
| `test_multi_partition_independent_computation`  | Two partitions with different F81 models compute independently |
| `test_multi_partition_internal_node_ab`         | Node AB posteriors differ between partitions                   |
| `test_multi_partition_sparse_dense_consistency` | Dense and sparse log-likelihoods match for clean sequences     |

---

## Consistency Tests

**File:** [`test_marginal_consistency.rs`](../../packages/treetime/src/commands/ancestral/__tests__/test_marginal_consistency.rs)

| Test                                                                     | Purpose                                                                 |
| ------------------------------------------------------------------------ | ----------------------------------------------------------------------- |
| `test_marginal_dense_sparse_log_lh_consistency_gap_free`                 | Dense and sparse log-likelihoods match (gap-free, JC69)                 |
| `test_marginal_sparse_varpos_matches_dense_profile_gap_free`             | Sparse variable-position distributions match dense rows                 |
| `test_marginal_dense_sparse_ambiguous_character_expectations_documented` | Log-likelihood agreement with ambiguity codes                           |
| `test_marginal_dense_sparse_ambiguous_r_reference_state_consistency`     | Dense and sparse reconstruct the same exact state for partial ambiguity |
| `test_marginal_posteriors_sum_to_one_skewed_gtr`                         | Normalization under extreme GTR (pi=[0.9,0.06,0.02,0.02])               |

---

## Branch Length Tests

**Files:** [`test_marginal_branch_length_equilibrium.rs`](../../packages/treetime/src/commands/ancestral/__tests__/test_marginal_branch_length/test_marginal_branch_length_equilibrium.rs), [`test_marginal_branch_length_monotonicity.rs`](../../packages/treetime/src/commands/ancestral/__tests__/test_marginal_branch_length/test_marginal_branch_length_monotonicity.rs)

### Equilibrium Convergence (7 tests)

| Test                                                  | Purpose                                                    |
| ----------------------------------------------------- | ---------------------------------------------------------- |
| `test_equilibrium_convergence_dense`                  | At t=100, log-likelihood converges to equilibrium product  |
| `test_equilibrium_convergence_sparse`                 | Same for sparse                                            |
| `test_equilibrium_convergence_nonuniform_pi`          | Equilibrium with pi=[0.4,0.1,0.2,0.3]                      |
| `test_equilibrium_convergence_multiple_positions`     | Multi-position equilibrium under JC69                      |
| `test_equilibrium_star_tree`                          | 4-leaf polytomy at equilibrium                             |
| `test_branch_length_sensitivity_near_zero`            | Numerical stability at t=1e-8                              |
| `test_dense_sparse_consistency_across_branch_lengths` | Dense/sparse match across [0.01, 0.1, 0.5, 1.0, 5.0, 20.0] |

### Monotonicity (5 tests)

| Test                                                        | Purpose                                                       |
| ----------------------------------------------------------- | ------------------------------------------------------------- |
| `test_likelihood_monotonic_increase_mismatched_sequences`   | Likelihood increases with branch length for mismatched leaves |
| `test_likelihood_maximized_near_zero_for_matched_sequences` | Likelihood decreases with branch length for matched leaves    |
| `test_likelihood_finite_across_branch_length_range`         | Log-likelihood finite and <=0 across [0.001, 10.0]            |
| `test_likelihood_monotonicity_three_taxon_tree`             | Monotonic decrease for identical sequences on 3-taxon tree    |
| `test_likelihood_monotonicity_sparse_partition`             | Sparse monotonicity for identical sequences                   |

---

## Topology Tests

**Files:** [`test_marginal_topology_caterpillar.rs`](../../packages/treetime/src/commands/ancestral/__tests__/test_marginal_topology/test_marginal_topology_caterpillar.rs), [`test_marginal_topology_deep_tree.rs`](../../packages/treetime/src/commands/ancestral/__tests__/test_marginal_topology/test_marginal_topology_deep_tree.rs), [`test_marginal_topology_polytomy.rs`](../../packages/treetime/src/commands/ancestral/__tests__/test_marginal_topology/test_marginal_topology_polytomy.rs)

### Caterpillar Trees (5 tests)

| Test                                             | Purpose                                                    |
| ------------------------------------------------ | ---------------------------------------------------------- |
| `test_caterpillar_tree_dense_completes`          | 5-leaf caterpillar produces finite log-likelihood          |
| `test_caterpillar_tree_sparse_completes`         | Same for sparse                                            |
| `test_caterpillar_tree_dense_sparse_consistency` | Dense and sparse match on varied sequences                 |
| `test_caterpillar_tree_varied_sequences`         | All-different sequences (AAAA, CCCC, GGGG, TTTT, ACGT)     |
| `test_caterpillar_tree_asymmetric_branches`      | Branch lengths [0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8] |

### Deep Trees (4 tests)

| Test                                                  | Purpose                                   |
| ----------------------------------------------------- | ----------------------------------------- |
| `test_deep_caterpillar_tree_10_leaves`                | 10-leaf deep tree (dense)                 |
| `test_deep_caterpillar_tree_sparse`                   | Same for sparse                           |
| `test_deep_caterpillar_tree_dense_sparse_consistency` | Dense/sparse match on 10-leaf tree        |
| `test_deep_tree_extreme_branches`                     | Extreme branch lengths [1e-8, 0.001, 5.0] |

### Polytomies (4 tests)

| Test                                     | Purpose                         |
| ---------------------------------------- | ------------------------------- |
| `test_polytomy_four_children`            | 4-way multifurcation at root    |
| `test_polytomy_dense_sparse_consistency` | Dense/sparse match for polytomy |
| `test_mixed_polytomy_binary`             | Mixed tree ((A,B,C)ABC,D,E)root |
| `test_large_polytomy`                    | 8-leaf polytomy                 |

---

## Stability Tests (Edge Cases)

**Files:** [`test_marginal_stability_extreme_branches.rs`](../../packages/treetime/src/commands/ancestral/__tests__/test_marginal_stability/test_marginal_stability_extreme_branches.rs), [`test_marginal_stability_near_zero_pi.rs`](../../packages/treetime/src/commands/ancestral/__tests__/test_marginal_stability/test_marginal_stability_near_zero_pi.rs), [`test_marginal_stability_rapid_transitions.rs`](../../packages/treetime/src/commands/ancestral/__tests__/test_marginal_stability/test_marginal_stability_rapid_transitions.rs)

### Extreme Branches (6 tests)

| Test                                      | Scenario                     |
| ----------------------------------------- | ---------------------------- |
| `test_extreme_short_branch_dense`         | t=1e-10, normalization holds |
| `test_extreme_short_branch_sparse`        | Same for sparse              |
| `test_extreme_long_branch_dense`          | t=10.0, normalization holds  |
| `test_extreme_long_branch_sparse`         | Same for sparse              |
| `test_extreme_asymmetric_branches_dense`  | t1=1e-10, t2=10.0            |
| `test_extreme_asymmetric_branches_sparse` | Same for sparse              |

### Near-Zero Pi (4 tests)

| Test                                          | Scenario                                           |
| --------------------------------------------- | -------------------------------------------------- |
| `test_near_zero_pi_dense`                     | pi=[0.97,0.01,0.01,0.01], observes rare states C/G |
| `test_near_zero_pi_sparse`                    | Same for sparse                                    |
| `test_near_zero_pi_with_dominant_state_dense` | Same pi, observes dominant state A                 |
| `test_extremely_skewed_pi_dense`              | pi=[0.9997,0.0001,0.0001,0.0001]                   |

### Rapid Transitions (6 tests)

| Test                                      | Scenario                                 |
| ----------------------------------------- | ---------------------------------------- |
| `test_high_mutation_rate_dense`           | mu=10.0                                  |
| `test_high_mutation_rate_sparse`          | Same for sparse                          |
| `test_very_high_mutation_rate_dense`      | mu=100.0                                 |
| `test_high_mutation_nonuniform_pi_dense`  | mu=10.0, pi=[0.4,0.1,0.2,0.3]            |
| `test_combined_extreme_parameters_dense`  | mu=50.0, pi=[0.9,0.03,0.04,0.03], t=1e-8 |
| `test_combined_extreme_parameters_sparse` | Same for sparse                          |

---

## Analytical Verification Tests

**Files:** [`test_marginal_analytical_two_taxon.rs`](../../packages/treetime/src/commands/ancestral/__tests__/test_marginal_analytical/test_marginal_analytical_two_taxon.rs), [`test_marginal_analytical_three_taxon.rs`](../../packages/treetime/src/commands/ancestral/__tests__/test_marginal_analytical/test_marginal_analytical_three_taxon.rs), [`test_marginal_analytical_star_tree.rs`](../../packages/treetime/src/commands/ancestral/__tests__/test_marginal_analytical/test_marginal_analytical_star_tree.rs)

### Two-Taxon (5 tests)

| Test                                              | Scenario                          |
| ------------------------------------------------- | --------------------------------- |
| `test_two_taxon_analytical_jc69_same_state`       | Both leaves A, analytical formula |
| `test_two_taxon_analytical_jc69_different_states` | Leaves A and T                    |
| `test_two_taxon_analytical_nonuniform_pi`         | GTR pi=[0.4,0.1,0.2,0.3]          |
| `test_two_taxon_analytical_multiple_positions`    | 3-position alignment              |
| `test_two_taxon_analytical_asymmetric_branches`   | t1=0.01, t2=1.0                   |

### Three-Taxon (2 tests)

| Test                                          | Scenario                                                 |
| --------------------------------------------- | -------------------------------------------------------- |
| `test_three_taxon_single_position_exhaustive` | Manual summation over all 16 internal state combinations |
| `test_three_taxon_all_combinations`           | All 64 single-position state triplets verified           |

### Star Tree (3 tests)

| Test                                          | Scenario                 |
| --------------------------------------------- | ------------------------ |
| `test_star_tree_analytical_jc69_all_same`     | 4 leaves all A           |
| `test_star_tree_analytical_jc69_mixed_states` | Leaves A, C, G, T        |
| `test_star_tree_analytical_nonuniform_pi`     | GTR pi=[0.1,0.2,0.3,0.4] |

**Purpose:** Verify implementation against closed-form Felsenstein likelihood formulas

---

## Softmax with Log-Norm Tests

**File:** [`softmax_with_log_norm.rs`](../../packages/treetime-utils/src/array/softmax_with_log_norm.rs) (inline `#[cfg(test)]`)

| Test                                         | Purpose                                                |
| -------------------------------------------- | ------------------------------------------------------ |
| `test_softmax_with_log_norm_finite`          | All-finite inputs: softmax has no zeros, all positive  |
| `test_softmax_with_log_norm_with_neg_inf`    | Mixed finite/-inf: only finite states get probability  |
| `test_softmax_with_log_norm_degenerate`      | All-NEG_INFINITY fallback: uniform, log_norm=-inf      |
| `test_softmax_with_log_norm_uniform_input`   | Equal finite inputs: exact uniform output              |
| `test_softmax_with_log_norm_shift_invariant` | Shift invariance: constant offset preserves softmax    |
| `stress::test_softmax_with_log_norm_*`       | Numerical stability: overflow/underflow boundary tests |

**Algorithm:** Fused logsumexp + softmax for numerically stable normalization. Subtracts max before exponentiation to prevent overflow/underflow. Returns both softmax result and log normalization constant.

---

## Dense Normalize-from-Log Tests

**File:** [`marginal_dense.rs`](../../packages/treetime/src/representation/partition/marginal_dense.rs) (inline `#[cfg(test)]`)

| Test                                                      | Purpose                                                            |
| --------------------------------------------------------- | ------------------------------------------------------------------ |
| `test_normalize_from_log_equal_probs`                     | Equal log-probs produce uniform distribution, finite log_lh        |
| `test_normalize_from_log_descending`                      | Unequal log-probs: analytically verified output and norm           |
| `test_normalize_from_log_all_neg_inf_single_row`          | All-NEG_INFINITY single row: uniform 1/4, log_lh=-inf              |
| `test_normalize_from_log_all_neg_inf_multiple_rows`       | All-NEG_INFINITY across 3 rows: uniform fallback each              |
| `test_normalize_from_log_mixed_neg_inf_and_finite_rows`   | Degenerate row 0 + normal row 1: row-level fallback, log_lh=-inf   |
| `test_normalize_from_log_mixed_finite_neg_inf_within_row` | Mixed finite/-inf within a row: only finite states get probability |
| `test_normalize_from_log_large_negative_values`           | Numerical stability: offset -1000 matches reference                |
| `test_normalize_from_log_three_states`                    | All-NEG_INFINITY fallback with 3 states: uniform 1/3               |
| `test_normalize_inplace_normal_rows`                      | Normal rows: correct normalization and log-likelihood              |
| `test_normalize_inplace_zero_row_returns_uniform`         | All-zero row: uniform fallback, log_lh=-inf                        |
| `test_normalize_inplace_mixed_zero_and_normal_rows`       | Zero row + normal row: per-row handling, log_lh=-inf               |
| `test_normalize_inplace_nan_row_returns_uniform`          | NaN row: uniform fallback, log_lh=-inf                             |
| `test_normalize_inplace_inf_row_returns_uniform`          | Inf row: uniform fallback, log_lh=-inf                             |

**Algorithm:** 2D logsumexp normalization matching v0's `normalize_profile(log=True)`. Guards all-`-inf` rows with uniform fallback, consistent with sparse path's `logsumexp_normalize`. `normalize_inplace` guards zero-sum and non-finite rows with the same uniform fallback.

---

## Substitution Composition Tests

**File:** [`test_mutation.rs`](../../packages/treetime/src/seq/__tests__/test_mutation.rs)

| Test                                                            | Purpose                                          |
| --------------------------------------------------------------- | ------------------------------------------------ |
| `test_mutation_compose_substitutions_both_empty`                | Empty parent and child produce empty result      |
| `test_mutation_compose_substitutions_parent_empty`              | Empty parent passes child through                |
| `test_mutation_compose_substitutions_child_empty`               | Empty child passes parent through                |
| `test_mutation_compose_substitutions_non_overlapping`           | Disjoint positions merge sorted                  |
| `test_mutation_compose_substitutions_chain`                     | A->G + G->T = A->T at same position              |
| `test_mutation_compose_substitutions_cancellation`              | A->G + G->A = no mutation                        |
| `test_mutation_compose_substitutions_mixed`                     | Chain, passthrough, and cancellation in one call |
| `test_mutation_compose_substitutions_output_sorted_by_position` | Interleaved positions verify merge order         |
| `test_mutation_compose_substitutions_all_cancel`                | All positions cancel, empty result               |

**Algorithm:** Two-pointer merge composition of parent and child substitution lists by position: chain, cancellation, or passthrough for non-overlapping positions.

### Prune integration tests for composition

**File:** [`test_run.rs`](../../packages/treetime/src/commands/prune/__tests__/test_run.rs)

| Test                                             | Purpose                                                      |
| ------------------------------------------------ | ------------------------------------------------------------ |
| `test_collapse_edge_compose_non_overlapping`     | Non-overlapping subs preserved through edge collapse         |
| `test_collapse_edge_compose_chain`               | Chain composition (A->G + G->T = A->T) through edge collapse |
| `test_collapse_edge_compose_cancellation`        | Cancellation (A->G + G->A = none) through edge collapse      |
| `test_collapse_edge_compose_multiple_partitions` | Composition applied independently per partition              |

---

## Property Test Generator Validation

**Directory:** [`prop_generators/`](../../packages/treetime/src/commands/ancestral/__tests__/prop_generators/)

### alignment.rs (4 tests)

| Test                                              | Purpose                                          |
| ------------------------------------------------- | ------------------------------------------------ |
| `test_prop_alignment_arb_sequence_length`         | Generated sequences have correct length          |
| `test_prop_alignment_arb_sequence_valid_chars`    | Characters are valid nucleotides/ambiguity codes |
| `test_prop_alignment_arb_alignment_structure`     | Alignment has correct taxa names and count       |
| `test_prop_alignment_arb_alignment_no_gaps_chars` | Gap-free alignment contains only ACGT            |

### branch_length.rs (2 tests)

| Test                                                 | Purpose                            |
| ---------------------------------------------------- | ---------------------------------- |
| `test_prop_branch_length_arb_branch_length_positive` | Branch lengths >= 0.001 and finite |
| `test_prop_branch_length_arb_branch_length_bounded`  | Branch lengths <= 2.0              |

### input.rs (3 tests)

| Test                                                          | Purpose                                          |
| ------------------------------------------------------------- | ------------------------------------------------ |
| `test_prop_input_arb_marginal_input_valid`                    | Newick, alignment size, sequence lengths, pi sum |
| `test_prop_input_arb_marginal_input_taxa_match`               | All alignment taxa appear in Newick              |
| `test_prop_input_arb_marginal_input_parseable_and_taxa_exact` | Parsed tree leaves match alignment taxa exactly  |

### tree.rs (6 tests)

| Test                                                       | Purpose                                       |
| ---------------------------------------------------------- | --------------------------------------------- |
| `test_prop_tree_arb_newick_valid_syntax`                   | Semicolon, taxa present, balanced parentheses |
| `test_prop_tree_arb_newick_has_branch_lengths`             | Newick contains colon branch lengths          |
| `test_prop_tree_arb_newick_all_taxa_present`               | All 5 taxa appear exactly once                |
| `test_prop_tree_arb_newick_no_double_branch_lengths`       | No ":0.1:0.2" patterns                        |
| `test_prop_tree_arb_newick_no_single_element_parens`       | No "(X)" without comma inside                 |
| `test_prop_tree_arb_newick_parseable_and_leaf_names_exact` | Parsed graph has correct leaf names           |

---

## Support Files (No Tests)

| File                                                                                                                                                           | Purpose                                                                              |
| -------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------ |
| [`prop_marginal_support.rs`](../../packages/treetime/src/commands/ancestral/__tests__/prop_marginal_support.rs)                                                | `run_dense_marginal()`, `run_sparse_marginal()` for property tests                   |
| [`test_marginal_analytical_support.rs`](../../packages/treetime/src/commands/ancestral/__tests__/test_marginal_analytical/test_marginal_analytical_support.rs) | Analytical likelihood formulas and `run_dense_marginal_get_log_lh()`                 |
| [`test_marginal_stability_support.rs`](../../packages/treetime/src/commands/ancestral/__tests__/test_marginal_stability/test_marginal_stability_support.rs)    | `assert_dense_profile_stable()`, `assert_sparse_profile_stable()`, partition runners |

---

## Test Tolerances

| Category                 | Tolerance       | Rationale                              |
| ------------------------ | --------------- | -------------------------------------- |
| Normalization            | 4 ULPs          | Machine precision for probability sums |
| Log-likelihood           | 1e-10 epsilon   | Accumulates floating-point error       |
| Dense/sparse             | 7.6e-6 relative | Fitch ambiguity resolution differences |
| Root invariance (sparse) | 1e-1            | Compression pattern varies with root   |
| Python parity            | 1e-6            | BLAS differences between NumPy/ndarray |
| Analytical               | 1e-7            | Closed-form formula comparison         |
