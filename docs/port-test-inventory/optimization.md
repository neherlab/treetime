# Branch Length Optimization Tests

[Back to index](_index.md)

## Summary

| Category                        | Files  | Tests  | Support Files | Type |
| ------------------------------- | ------ | ------ | ------------- | ---- |
| Coefficient extraction (dense)  | 5      | 13     | 1             | Unit |
| Coefficient extraction (sparse) | 7      | 19     | 0             | Unit |
| Newton-Raphson convergence      | 2      | 5      | 1             | Unit |
| Grid search                     | 3      | 8      | 1             | Unit |
| Dense/sparse equivalence        | 4      | 8      | 1             | Unit |
| Convergence control             | 3      | 8      | 1             | Unit |
| Optimization metrics            | 1      | 7      | 0             | Unit |
| Zero branch optimal             | 1      | 13     | 0             | Unit |
| Initial guess GTR messages      | 1      | 2      | 0             | Unit |
| **Total**                       | **27** | **83** | **5**         | Unit |

---

## Coefficient Extraction - Dense

**Directory:** [`test_coefficient_extraction_dense/`](../../packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_dense/)

### Basic Tests

**File:** [`test_coefficient_extraction_dense_basic.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_dense/test_coefficient_extraction_dense_basic.rs)

| Test                                              | Purpose                          |
| ------------------------------------------------- | -------------------------------- |
| `test_get_coefficients_identity_messages`         | Uniform messages coefficient sum |
| `test_get_coefficients_certain_state_parent`      | Certain parent, uniform child    |
| `test_get_coefficients_certain_state_child`       | Uniform parent, certain child    |
| `test_get_coefficients_matching_certain_states`   | Both certain same state          |
| `test_get_coefficients_mismatched_certain_states` | Certain A vs certain C mismatch  |

### GTR Model Tests

**File:** [`test_coefficient_extraction_dense_gtr.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_dense/test_coefficient_extraction_dense_gtr.rs)

| Test                                              | Purpose                          |
| ------------------------------------------------- | -------------------------------- |
| `test_coefficients_use_eigenvector_decomposition` | Verify eigenvector-based formula |
| `test_coefficients_preserve_gtr_reference`        | Stored GTR eigenvalues match     |

### Partition Tests

**File:** [`test_coefficient_extraction_dense_partition.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_dense/test_coefficient_extraction_dense_partition.rs)

| Test                              | Purpose                           |
| --------------------------------- | --------------------------------- |
| `test_partition_contribution_new` | PartitionContribution stores data |

### Position Tests

**File:** [`test_coefficient_extraction_dense_positions.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_dense/test_coefficient_extraction_dense_positions.rs)

| Test                                       | Purpose                          |
| ------------------------------------------ | -------------------------------- |
| `test_get_coefficients_multiple_positions` | Multi-position coefficient shape |
| `test_get_coefficients_row_independence`   | Positions computed independently |

### Properties Tests

**File:** [`test_coefficient_extraction_dense_properties.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_dense/test_coefficient_extraction_dense_properties.rs)

| Test                                                          | Purpose                            |
| ------------------------------------------------------------- | ---------------------------------- |
| `test_coefficients_produce_valid_likelihood_at_zero`          | Finite log-LH at zero branch       |
| `test_coefficients_likelihood_decreases_for_mismatch_at_zero` | Mismatch lower LH than match       |
| `test_coefficients_support_positive_branch_length`            | Valid metrics at positive branches |

### Support

**File:** [`test_coefficient_extraction_dense_support.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_dense/test_coefficient_extraction_dense_support.rs)

Helper file. Provides `make_dense_seq_dis()`. No tests.

---

## Coefficient Extraction - Sparse

**Directory:** [`test_coefficient_extraction_sparse/`](../../packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_sparse/)

### Site Tests

**File:** [`test_coefficient_extraction_sparse_site.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_sparse/test_coefficient_extraction_sparse_site.rs)

| Test                                                          | Purpose                         |
| ------------------------------------------------------------- | ------------------------------- |
| `test_site_contribution_stores_multiplicity_and_coefficients` | Store multiplicity and coeffs   |
| `test_site_contribution_unit_multiplicity`                    | Unit multiplicity for variables |

### Coefficients Tests

**File:** [`test_coefficient_extraction_sparse_coefficients.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_sparse/test_coefficient_extraction_sparse_coefficients.rs)

| Test                                                       | Purpose                      |
| ---------------------------------------------------------- | ---------------------------- |
| `test_coefficients_computed_via_eigenvector_decomposition` | Eigenvector-based formula    |
| `test_matching_states_high_lh_at_zero`                     | High LH for matching states  |
| `test_mismatched_states_low_lh_at_zero`                    | Low LH for mismatched states |

### Derivatives Tests

**File:** [`test_coefficient_extraction_sparse_derivatives.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_sparse/test_coefficient_extraction_sparse_derivatives.rs)

| Test                                                       | Purpose                             |
| ---------------------------------------------------------- | ----------------------------------- |
| `test_evaluate_sparse_derivative_matches_numerical`        | Analytical vs numerical derivative  |
| `test_evaluate_sparse_derivative_scales_with_multiplicity` | Derivative scales with multiplicity |

### Eigenvalues Tests

**File:** [`test_coefficient_extraction_sparse_eigenvalues.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_sparse/test_coefficient_extraction_sparse_eigenvalues.rs)

| Test                                               | Purpose                          |
| -------------------------------------------------- | -------------------------------- |
| `test_eigenvalues_affect_branch_length_evaluation` | LH differs at different branches |
| `test_jc69_eigenvalues_structure`                  | One zero, three negative eigvals |

### Evaluate Tests

**File:** [`test_coefficient_extraction_sparse_evaluate.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_sparse/test_coefficient_extraction_sparse_evaluate.rs)

| Test                                               | Purpose                     |
| -------------------------------------------------- | --------------------------- |
| `test_evaluate_sparse_single_site_multiplicity_1`  | Single site LH at zero      |
| `test_evaluate_sparse_single_site_multiplicity_10` | Multiplicity 10 LH at zero  |
| `test_evaluate_sparse_multiplicity_scales_log_lh`  | LH scales with multiplicity |
| `test_evaluate_sparse_multiple_sites_sum_log_lh`   | Multiple sites sum log-LH   |

### Mixed Sites Tests

**File:** [`test_coefficient_extraction_sparse_mixed_sites.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_sparse/test_coefficient_extraction_sparse_mixed_sites.rs)

| Test                                               | Purpose                        |
| -------------------------------------------------- | ------------------------------ |
| `test_mixed_variable_and_fixed_sites`              | Variable + fixed sites compose |
| `test_fixed_sites_dominate_with_high_multiplicity` | High multiplicity dominates    |

### Partition Tests

**File:** [`test_coefficient_extraction_sparse_partition.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_sparse/test_coefficient_extraction_sparse_partition.rs)

| Test                                                        | Purpose                       |
| ----------------------------------------------------------- | ----------------------------- |
| `test_partition_contribution_empty_sites`                   | Empty site contributions      |
| `test_partition_contribution_single_variable_site`          | Single variable site          |
| `test_partition_contribution_multiple_variable_sites`       | Multiple variable sites       |
| `test_partition_contribution_fixed_sites_with_multiplicity` | Fixed sites with multiplicity |

---

## Newton-Raphson Convergence

**Directory:** [`test_newton_convergence/`](../../packages/treetime/src/commands/optimize/__tests__/test_newton_convergence/)

### Iteration Tests

**File:** [`test_newton_convergence_iteration.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_newton_convergence/test_newton_convergence_iteration.rs)

| Test                                            | Purpose                   |
| ----------------------------------------------- | ------------------------- |
| `test_newton_iteration_converges_within_bounds` | Convergence within bounds |
| `test_newton_iteration_respects_max_iter`       | Stops at max iterations   |

### Evaluate Mixed Tests

**File:** [`test_newton_convergence_evaluate_mixed.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_newton_convergence/test_newton_convergence_evaluate_mixed.rs)

| Test                                          | Purpose                       |
| --------------------------------------------- | ----------------------------- |
| `test_evaluate_mixed_returns_finite_values`   | Finite values for valid input |
| `test_evaluate_mixed_log_lh_negative`         | Log-LH is non-positive        |
| `test_evaluate_mixed_multiple_branch_lengths` | Finite across branch lengths  |

### Support

**File:** [`test_newton_convergence_support.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_newton_convergence/test_newton_convergence_support.rs)

Helper file. Provides `make_dense_contribution()`. No tests.

---

## Grid Search

**Directory:** [`test_grid_search/`](../../packages/treetime/src/commands/optimize/__tests__/test_grid_search/)

### Basic Tests

**File:** [`test_grid_search_basic.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_grid_search/test_grid_search_basic.rs)

| Test                                        | Purpose                     |
| ------------------------------------------- | --------------------------- |
| `test_grid_search_finds_maximum_log_lh`     | Finds likelihood maximum    |
| `test_grid_search_result_in_range`          | Result within search bounds |
| `test_grid_search_with_small_branch_length` | Works with small branches   |
| `test_grid_search_with_large_branch_length` | Works with large branches   |

### Edge Cases Tests

**File:** [`test_grid_search_edge_cases.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_grid_search/test_grid_search_edge_cases.rs)

| Test                                                | Purpose                      |
| --------------------------------------------------- | ---------------------------- |
| `test_grid_search_one_mutation_dominates_range`     | Large one_mutation vs branch |
| `test_grid_search_preserves_positive_branch_length` | Always returns positive      |

### Partitions Tests

**File:** [`test_grid_search_partitions.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_grid_search/test_grid_search_partitions.rs)

| Test                                    | Purpose                         |
| --------------------------------------- | ------------------------------- |
| `test_grid_search_multiple_partitions`  | Combines multiple contributions |
| `test_grid_search_uniform_coefficients` | Handles uniform coefficients    |

### Support

**File:** [`test_grid_search_support.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_grid_search/test_grid_search_support.rs)

Helper file. Provides `make_dense_contribution()` and `grid_search()`. No tests.

---

## Dense/Sparse Equivalence

**Directory:** [`test_dense_sparse_equivalence/`](../../packages/treetime/src/commands/optimize/__tests__/test_dense_sparse_equivalence/)

### Initial Tests

**File:** [`test_dense_sparse_equivalence_initial.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_initial.rs)

| Test                                                          | Purpose                      |
| ------------------------------------------------------------- | ---------------------------- |
| `test_dense_sparse_initial_log_lh_equivalence`                | Initial LH matches           |
| `test_dense_sparse_initial_log_lh_equivalence_with_mutations` | Initial LH matches with muts |

### Convergence Tests

**File:** [`test_dense_sparse_equivalence_convergence.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_convergence.rs)

| Test                                 | Purpose               |
| ------------------------------------ | --------------------- |
| `test_dense_optimization_converges`  | Dense mode converges  |
| `test_sparse_optimization_converges` | Sparse mode converges |

### Bounds Tests

**File:** [`test_dense_sparse_equivalence_bounds.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_bounds.rs)

| Test                                                             | Purpose                    |
| ---------------------------------------------------------------- | -------------------------- |
| `test_dense_sparse_log_lh_bounded_difference_after_optimization` | LH difference is bounded   |
| `test_dense_sparse_branch_lengths_bounded_difference`            | Branch lengths are similar |

### Validity Tests

**File:** [`test_dense_sparse_equivalence_validity.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_validity.rs)

| Test                                              | Purpose                       |
| ------------------------------------------------- | ----------------------------- |
| `test_dense_optimization_produces_valid_results`  | Dense produces valid results  |
| `test_sparse_optimization_produces_valid_results` | Sparse produces valid results |

### Support

**File:** [`test_dense_sparse_equivalence_support.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_support.rs)

Helper file. Provides `setup_dense_only()`, `setup_sparse_only()`, `get_branch_lengths()`, `gap_free_alignment()`. No tests.

---

## Convergence Control

**Directory:** [`test_convergence/`](../../packages/treetime/src/commands/optimize/__tests__/test_convergence/)

### Iterations Tests

**File:** [`test_convergence_iterations.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_convergence/test_convergence_iterations.rs)

| Test                                                 | Purpose                          |
| ---------------------------------------------------- | -------------------------------- |
| `test_optimization_converges_within_iterations`      | Converges within iteration limit |
| `test_optimization_improves_or_maintains_likelihood` | LH improves or stays stable      |
| `test_optimization_produces_valid_branch_lengths`    | Branch lengths valid after opt   |

### Edge Cases Tests

**File:** [`test_convergence_edge_cases.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_convergence/test_convergence_edge_cases.rs)

| Test                                            | Purpose                      |
| ----------------------------------------------- | ---------------------------- |
| `test_optimization_handles_zero_branch_lengths` | Handles zero-length branches |
| `test_optimization_handles_very_short_branches` | Handles very short branches  |
| `test_optimization_handles_long_branches`       | Handles long branches        |

### Idempotence Tests

**File:** [`test_convergence_idempotence.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_convergence/test_convergence_idempotence.rs)

| Test                                                    | Purpose                        |
| ------------------------------------------------------- | ------------------------------ |
| `test_optimization_converges_with_valid_branch_lengths` | Branch lengths valid each iter |
| `test_second_optimization_produces_same_likelihood`     | Deterministic across runs      |

### Support

**File:** [`test_convergence_support.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_convergence/test_convergence_support.rs)

Helper file. Provides `simple_alignment()`, `setup_partitions()`, `compute_total_lh()`. No tests.

---

## Optimization Metrics

**File:** [`test_optimization_metrics.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_optimization_metrics.rs)

| Test                                                 | Purpose                       |
| ---------------------------------------------------- | ----------------------------- |
| `test_optimization_metrics_default_is_zero`          | Default metrics are zero      |
| `test_optimization_metrics_new`                      | Constructor stores values     |
| `test_optimization_metrics_add_single`               | Add single partition          |
| `test_optimization_metrics_add_multiple_partitions`  | Add multiple partitions       |
| `test_optimization_metrics_add_with_negative_values` | Add with negative values      |
| `test_optimization_metrics_add_preserves_precision`  | Precision with 1000 additions |
| `test_optimization_metrics_add_zero_is_identity`     | Adding zero is identity       |

---

## Zero Branch Optimal

**File:** [`test_is_zero_branch_optimal.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_is_zero_branch_optimal.rs)

| Test                                                                    | Purpose                          |
| ----------------------------------------------------------------------- | -------------------------------- |
| `test_is_zero_branch_optimal_empty_contributions`                       | Empty contributions return false |
| `test_is_zero_branch_optimal_high_lh_negative_derivative_returns_true`  | High LH + negative deriv = true  |
| `test_is_zero_branch_optimal_high_lh_positive_derivative_returns_false` | High LH + positive deriv = false |
| `test_is_zero_branch_optimal_low_lh_skips_derivative`                   | Low LH skips derivative check    |
| `test_is_zero_branch_optimal_multiple_partitions_product_lh`            | Multiple partitions product LH   |
| `test_is_zero_branch_optimal_multiple_partitions_low_combined_lh`       | Low combined LH returns false    |
| `test_is_zero_branch_optimal_sparse_high_lh_negative_derivative`        | Sparse with negative derivative  |
| `test_is_zero_branch_optimal_sparse_high_lh_zero_derivative`            | Sparse with zero derivative      |
| `test_is_zero_branch_optimal_mixed_dense_and_sparse`                    | Mixed dense and sparse           |
| `test_is_zero_branch_optimal_multiple_positions_dense`                  | Multiple dense positions         |
| `test_is_zero_branch_optimal_sparse_multiplicity_effect`                | Multiplicity affects result      |
| `test_is_zero_branch_optimal_boundary_lh_just_above_threshold`          | LH just above 0.01 threshold     |
| `test_is_zero_branch_optimal_boundary_lh_at_threshold`                  | LH at exact threshold boundary   |

---

## Initial Guess GTR Messages

**File:** [`test_initial_guess_gtr_messages.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_initial_guess_gtr_messages.rs)

Regression tests verifying that `initial_guess_mixed()` reads edge messages computed with the real GTR model, not stale JC69 messages from the dummy initialization pass.

| Test                                             | Purpose                                                       |
| ------------------------------------------------ | ------------------------------------------------------------- |
| `test_stale_jc69_messages_bias_initial_guess`    | Stale JC69 messages produce different branch lengths than F81 |
| `test_initial_guess_idempotent_after_gtr_update` | Identical initialization sequences produce identical results  |
