# Branch Length Optimization Tests

[Back to index](_index.md)

## Summary

| Category                         | Files  | Tests   | Support Files | Type     |
| -------------------------------- | ------ | ------- | ------------- | -------- |
| Coefficient extraction (dense)   | 5      | 13      | 1             | Unit     |
| Coefficient extraction (sparse)  | 7      | 19      | 0             | Unit     |
| Newton-Raphson convergence       | 2      | 5       | 1             | Unit     |
| Grid search                      | 3      | 8       | 1             | Unit     |
| Dense/sparse equivalence         | 4      | 8       | 1             | Unit     |
| Convergence control              | 3      | 8       | 1             | Unit     |
| Optimization metrics             | 1      | 7       | 0             | Unit     |
| Zero branch optimal              | 1      | 13      | 0             | Unit     |
| Initial guess GTR messages       | 1      | 2       | 0             | Unit     |
| Initial guess soft Hamming       | 1      | 10      | 0             | Unit     |
| Topology cleanup in loop         | 1      | 15      | 0             | Unit     |
| Indel contribution (inline)      | 1      | 8       | 0             | Unit     |
| Indel contribution (integration) | 1      | 11      | 0             | Unit     |
| Indel contribution (property)    | 1      | 3       | 0             | Property |
| **Total**                        | **32** | **130** | **5**         |          |

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

| Test                                                                 | Purpose                                         |
| -------------------------------------------------------------------- | ----------------------------------------------- |
| `test_evaluate_sparse_first_derivative_matches_numerical`            | Analytical vs numerical first derivative        |
| `test_evaluate_sparse_second_derivative_matches_numerical`           | Analytical vs numerical second derivative       |
| `test_evaluate_sparse_second_derivative_numerical_high_multiplicity` | Second derivative numerical at multiplicity=500 |
| `test_evaluate_sparse_derivatives_scale_with_multiplicity`           | All metrics scale linearly with multiplicity    |
| `test_evaluate_sparse_matches_dense_multiplicity_1`                  | Sparse (m=1) matches dense for all metrics      |

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

### Tolerance Tests

**File:** [`test_newton_convergence_tolerance.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_newton_convergence/test_newton_convergence_tolerance.rs)

| Test                                                     | Purpose                                    |
| -------------------------------------------------------- | ------------------------------------------ |
| `test_newton_tolerance_positive_branch_uses_relative`    | Relative tolerance dominates for large BL  |
| `test_newton_tolerance_zero_branch_uses_absolute_floor`  | Absolute floor prevents degenerate zero    |
| `test_newton_tolerance_small_branch_uses_absolute_floor` | Absolute floor for very small BL           |
| `test_newton_tolerance_crossover`                        | Relative = absolute at crossover point     |
| `test_newton_tolerance_always_positive`                  | Positive and finite for all BL values      |
| `test_newton_zero_start_converges_efficiently`           | Zero-start does not exhaust all iterations |

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

### Coverage Tests

**File:** [`test_grid_search_coverage.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_grid_search/test_grid_search_coverage.rs)

| Test                                                                  | Purpose                                                  |
| --------------------------------------------------------------------- | -------------------------------------------------------- |
| `test_grid_search_branch_lengths_zero_branch_covers_full_range`       | Grid extends to 0.5 subs/site when branch_length=0       |
| `test_grid_search_branch_lengths_large_branch_extends_beyond_minimum` | Proportional bound dominates for large branches          |
| `test_grid_search_branch_lengths_is_log_spaced`                       | Grid is log-spaced with constant ratio                   |
| `test_grid_search_branch_lengths_has_correct_point_count`             | Grid always has 100 points                               |
| `test_grid_search_reaches_extended_range_from_zero`                   | Monotonic likelihood finds optimum near upper bound      |
| `test_grid_search_branch_lengths_coverage_invariant`                  | Grid covers [one_mutation, 0.5] for all parameter combos |

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

| Test                                                 | Purpose                                                   |
| ---------------------------------------------------- | --------------------------------------------------------- |
| `test_optimization_converges_within_iterations`      | Bounded oscillation: undamped LH stays within tight range |
| `test_optimization_improves_or_maintains_likelihood` | Strict non-regression: final LH >= initial LH             |
| `test_optimization_produces_valid_branch_lengths`    | Branch lengths valid after opt                            |

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

## Outer-Loop Damping

**File:** [`test_damping.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_damping.rs)

| Test                                                     | Cases | Purpose                                                                           |
| -------------------------------------------------------- | ----- | --------------------------------------------------------------------------------- |
| `test_save_branch_lengths_captures_all_edges`            | 1     | Verifies save_branch_lengths reads all edges                                      |
| `test_apply_damping_zero_is_noop`                        | 1     | damping=0.0 leaves optimized values unchanged                                     |
| `test_apply_damping_weights_match_v0`                    | 5     | Blend weights match v0 formula at iterations 0-9                                  |
| `test_apply_damping_blends_correctly`                    | 1     | Verifies bl = new*w_new + old*w_old                                               |
| `test_apply_damping_new_weight_increases_with_iteration` | 1     | New weight monotonically increases over iterations                                |
| `test_damped_optimization_converges`                     | 1     | Production convergence criterion, sign-flip oscillation detection, tail stability |
| `test_damped_optimization_does_not_regress`              | 1     | Strict non-regression: final LH >= initial LH                                     |

---

## v0 Parity and Damped vs Undamped (Golden Master)

**File:** [`test_gm_optimize.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_gm_optimize.rs)

| Test                                  | Cases | Purpose                                                        |
| ------------------------------------- | ----- | -------------------------------------------------------------- |
| `test_gm_optimize`                    | 1     | v0/v1 total branch length parity within 5% (flu/h3n2/20, JC69) |
| `test_gm_optimize_damped_vs_undamped` | 1     | Damping does not increase oscillation on real dataset          |

Fixtures in `__fixtures__/`:

- `gm_optimize_inputs.json` - shared test parameters
- `gm_optimize_outputs.json` - v0 reference captured by `gm_optimize_capture`
- `gm_optimize_capture` - Python capture script for v0 TreeAnc.optimize_tree_marginal

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

## Initial Guess Mode

**File:** [`test_initial_guess_mode.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_initial_guess_mode.rs)

Tests for `InitialGuessMode` dispatch: NaN detection, selective fill (auto), overwrite (always), and error on missing (never).

| Test                                                    | Purpose                                                      |
| ------------------------------------------------------- | ------------------------------------------------------------ |
| `test_initial_guess_mode_default_is_auto`               | Default variant is Auto                                      |
| `test_initial_guess_mode_detects_nan_from_newick`       | Detects NaN branch lengths from bio crate newick parser      |
| `test_initial_guess_mode_detects_explicit_nan`          | Detects explicit `Some(NaN)` as missing                      |
| `test_initial_guess_mode_no_missing_when_all_finite`    | Reports no missing when all edges have finite lengths        |
| `test_initial_guess_mode_auto_preserves_valid_lengths`  | Auto mode does not overwrite valid finite branch lengths     |
| `test_initial_guess_mode_auto_fills_nan_from_newick`    | Auto mode fills NaN edges with finite non-negative values    |
| `test_initial_guess_mode_auto_fills_only_missing_edges` | Auto mode fills one NaN edge, preserves all others exactly   |
| `test_initial_guess_mode_always_overwrites_all`         | Always mode overwrites existing branch lengths               |
| `test_initial_guess_mode_never_accepts_complete_tree`   | Never mode accepts tree with all finite branch lengths       |
| `test_initial_guess_mode_never_rejects_nan_tree`        | Never mode detects NaN tree as having missing branch lengths |

---

## Initial Guess GTR Messages

**File:** [`test_initial_guess_gtr_messages.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_initial_guess_gtr_messages.rs)

Regression tests verifying that `initial_guess_mixed()` reads edge messages computed with the real GTR model, not stale JC69 messages from the dummy initialization pass.

| Test                                             | Purpose                                                       |
| ------------------------------------------------ | ------------------------------------------------------------- |
| `test_stale_jc69_messages_bias_initial_guess`    | Stale JC69 messages produce different branch lengths than F81 |
| `test_initial_guess_idempotent_after_gtr_update` | Identical initialization sequences produce identical results  |

---

## Initial Guess Soft Hamming

**File:** [`test_initial_guess_soft_hamming.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_initial_guess_soft_hamming.rs)

Tests for the dense partition's soft Hamming distance used in `initial_guess_mixed()`. Validates that `edge_initial_differences()` computes `sum(1 - dot(pp, pc))` correctly across profile shapes.

### Direct Numerical Tests

| Test                                                 | Purpose                                            |
| ---------------------------------------------------- | -------------------------------------------------- |
| `test_sharp_matching_profiles_zero_differences`      | One-hot matching profiles: dot=1, diff=0           |
| `test_sharp_mismatched_profiles_integer_differences` | One-hot mismatched profiles: dot=0, diff=1 each    |
| `test_uniform_profiles_fractional_differences`       | Uniform 0.25: dot=0.25, diff=0.75 per position     |
| `test_weakly_informative_same_dominant_state`        | Dominant 0.7 same state: dot=0.52, diff=0.48       |
| `test_weakly_informative_different_dominant_state`   | Dominant 0.7 different states: dot=0.16, diff=0.84 |
| `test_mixed_profiles_sum_of_contributions`           | Sum of sharp, uniform, and weak contributions      |
| `test_gap_positions_excluded`                        | Gap ranges excluded from soft Hamming sum          |

### Pipeline Integration Tests

| Test                                              | Purpose                                            |
| ------------------------------------------------- | -------------------------------------------------- |
| `test_identical_sequences_hard_zero_soft_small`   | Identical sequences: hard=0, soft small positive   |
| `test_divergent_sequences_soft_differs_from_hard` | Divergent sequences: soft and hard disagree        |
| `test_differences_bounded_by_effective_length`    | 0 <= differences <= effective_length for all edges |

---

## Topology Cleanup in Loop

**File:** [`test_topology_cleanup.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_topology_cleanup.rs)

| Test                                                               | Purpose                                                    |
| ------------------------------------------------------------------ | ---------------------------------------------------------- |
| `test_optimize_find_zero_optimal_internal_edges_empty_graph`       | Empty graph returns no edges                               |
| `test_optimize_find_zero_optimal_internal_edges_no_zero_edges`     | Non-zero branches return no edges                          |
| `test_optimize_find_zero_optimal_internal_edges_skips_leaves`      | Leaf edges with bl=0 are not collected                     |
| `test_optimize_find_zero_optimal_internal_edges_collects_internal` | Internal edges with bl=0 are collected                     |
| `test_optimize_find_zero_optimal_internal_edges_multiple`          | Multiple zero-optimal internal edges collected             |
| `test_optimize_collapse_edge_sparse_composes_subs`                 | Substitutions composed correctly on sparse collapse        |
| `test_optimize_collapse_edge_dense_cleanup`                        | Stale dense partition data removed after collapse          |
| `test_optimize_collapse_edge_branch_length_sum`                    | Branch lengths summed correctly (0 + child = child)        |
| `test_optimize_prune_and_merge_empty_list`                         | Empty zero-edge list is a noop                             |
| `test_optimize_prune_and_merge_collapses_and_merges`               | Collapse + shared mutation merge in one pass               |
| `test_optimize_loop_with_topology_cleanup_sparse`                  | Full loop collapses zero-optimal branches (identical seqs) |
| `test_optimize_loop_no_collapse_when_branches_nonzero`             | No collapse when all branches carry signal                 |
| `test_optimize_merge_then_marginal_finite_likelihood`              | Merge + update_marginal produces finite log-likelihood     |
| `test_optimize_cascading_collapse_parent_child_both_zero`          | Parent-child both zero-optimal, guard handles removal      |
| `test_optimize_loop_with_topology_cleanup_dense`                   | Dense-mode full loop with topology cleanup                 |
