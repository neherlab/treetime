# Branch Length Optimization Tests

[Back to index](README.md)

## Summary

| Category                                                                                | Type            |
| --------------------------------------------------------------------------------------- | --------------- |
| [Coefficient extraction (dense)](#coefficient-extraction---dense)                       | Unit + Property |
| [Coefficient extraction (dense invariants)](#coefficient-extraction---dense-invariants) | Parameterized   |
| [Coefficient extraction (sparse)](#coefficient-extraction---sparse)                     | Unit            |
| [Newton-Raphson convergence](#newton-raphson-convergence)                               | Unit            |
| [Grid search](#grid-search)                                                             | Unit            |
| [Dense/sparse equivalence](#densesparse-equivalence)                                    | Unit            |
| [Convergence control](#convergence-control)                                             | Unit            |
| [Convergence conditions](#convergence-conditions)                                       | Unit            |
| [Convergence on real datasets](#convergence-on-real-datasets)                           | Unit            |
| [Optimization metrics](#optimization-metrics)                                           | Unit            |
| [Zero branch optimal](#zero-branch-optimal)                                             | Unit            |
| [Initial guess mode](#initial-guess-mode)                                               | Unit            |
| [Initial guess indel zero-BL](#initial-guess-indel-zero-bl)                             | Unit            |
| [Initial guess GTR messages](#initial-guess-gtr-messages)                               | Unit            |
| [Initial guess formula](#initial-guess-formula)                                         | Unit            |
| [Initial guess gaps](#initial-guess-gaps)                                               | Unit            |
| [Dense edge subs](#dense-edge-subs)                                                     | Unit            |
| [Eval zero-branch mismatch](#eval-zero-branch-mismatch)                                 | Unit            |
| [Topology cleanup in loop](#topology-cleanup-in-loop)                                   | Unit            |
| [Indel contribution (inline)](#indel-contribution-inline)                               | Unit            |
| [Indel contribution (integration)](#indel-contribution-integration)                     | Unit + Property |
| [Optimization method](#optimization-method)                                             | Unit            |
| [Optimization method step clamping](#optimization-method-step-clamping)                 | Unit            |
| [Dispatch zero boundary](#dispatch-zero-boundary)                                       | Unit            |
| [Outer-loop damping](#outer-loop-damping)                                               | Unit            |
| [v0 parity (golden master)](#v0-parity-and-damped-vs-undamped-golden-master)            | Golden-master   |
| [Run optimize loop contract](#run-optimize-loop-contract)                               | Unit            |
| [CLI args](#cli-args)                                                                   | Unit            |
| [Zero sequence length](#zero-sequence-length)                                           | Unit            |

---

## Coefficient Extraction - Dense

**Directory:** [`test_coefficient_extraction_dense/`](../../packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_dense/)

### Basic Tests

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_dense/test_coefficient_extraction_dense_basic.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_dense/test_coefficient_extraction_dense_basic.rs)

**Impl:**

- [`packages/treetime/src/commands/optimize/optimize_dense.rs`](../../packages/treetime/src/commands/optimize/optimize_dense.rs)
- [`packages/treetime/src/commands/optimize/optimize_dense_eval.rs`](../../packages/treetime/src/commands/optimize/optimize_dense_eval.rs)

| Test                                              | Purpose                          |
| ------------------------------------------------- | -------------------------------- |
| `test_get_coefficients_identity_messages`         | Uniform messages coefficient sum |
| `test_get_coefficients_certain_state_parent`      | Certain parent, uniform child    |
| `test_get_coefficients_certain_state_child`       | Uniform parent, certain child    |
| `test_get_coefficients_matching_certain_states`   | Both certain same state          |
| `test_get_coefficients_mismatched_certain_states` | Certain A vs certain C mismatch  |

### GTR Model Tests

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_dense/test_coefficient_extraction_dense_gtr.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_dense/test_coefficient_extraction_dense_gtr.rs)

**Impl:** [`packages/treetime/src/commands/optimize/optimize_dense.rs`](../../packages/treetime/src/commands/optimize/optimize_dense.rs)

| Test                                              | Purpose                          |
| ------------------------------------------------- | -------------------------------- |
| `test_coefficients_use_eigenvector_decomposition` | Verify eigenvector-based formula |
| `test_coefficients_preserve_gtr_reference`        | Stored GTR eigenvalues match     |

### Partition Tests

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_dense/test_coefficient_extraction_dense_partition.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_dense/test_coefficient_extraction_dense_partition.rs)

**Impl:** [`packages/treetime/src/commands/optimize/optimize_dense.rs`](../../packages/treetime/src/commands/optimize/optimize_dense.rs)

| Test                              | Purpose                           |
| --------------------------------- | --------------------------------- |
| `test_partition_contribution_new` | PartitionContribution stores data |

### Position Tests

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_dense/test_coefficient_extraction_dense_positions.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_dense/test_coefficient_extraction_dense_positions.rs)

**Impl:** [`packages/treetime/src/commands/optimize/optimize_dense.rs`](../../packages/treetime/src/commands/optimize/optimize_dense.rs)

| Test                                       | Purpose                          |
| ------------------------------------------ | -------------------------------- |
| `test_get_coefficients_multiple_positions` | Multi-position coefficient shape |
| `test_get_coefficients_row_independence`   | Positions computed independently |

### Properties Tests

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_dense/test_coefficient_extraction_dense_properties.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_dense/test_coefficient_extraction_dense_properties.rs)

**Impl:**

- [`packages/treetime/src/commands/optimize/optimize_dense.rs`](../../packages/treetime/src/commands/optimize/optimize_dense.rs)
- [`packages/treetime/src/commands/optimize/optimize_dense_eval.rs`](../../packages/treetime/src/commands/optimize/optimize_dense_eval.rs)

| Test                                                          | Purpose                            |
| ------------------------------------------------------------- | ---------------------------------- |
| `test_coefficients_produce_valid_likelihood_at_zero`          | Finite log-LH at zero branch       |
| `test_coefficients_likelihood_decreases_for_mismatch_at_zero` | Mismatch lower LH than match       |
| `test_coefficients_support_positive_branch_length`            | Valid metrics at positive branches |

### Property Invariants

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_dense/test_coefficient_extraction_dense_prop_invariants.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_dense/test_coefficient_extraction_dense_prop_invariants.rs)

**Impl:**

- [`packages/treetime/src/commands/optimize/optimize_dense.rs`](../../packages/treetime/src/commands/optimize/optimize_dense.rs)
- [`packages/treetime/src/commands/optimize/optimize_dense_eval.rs`](../../packages/treetime/src/commands/optimize/optimize_dense_eval.rs)
- [`packages/treetime/src/commands/optimize/optimize_sparse.rs`](../../packages/treetime/src/commands/optimize/optimize_sparse.rs)
- [`packages/treetime/src/commands/optimize/optimize_sparse_eval.rs`](../../packages/treetime/src/commands/optimize/optimize_sparse_eval.rs)

| Test                                                               | Purpose                                                                            | Notes                 |
| ------------------------------------------------------------------ | ---------------------------------------------------------------------------------- | --------------------- |
| `test_coefficient_boundary_disjoint_support_zero_at_t0`            | Disjoint-support singular boundary at t=0                                          |                       |
| `test_prop_coefficient_nonneg_site_lh_at_zero`                     | Site likelihood positive for overlapping probability vectors                       |                       |
| `test_prop_coefficient_multiplicity_linearity`                     | Sparse multiplicity factor acts linearly on all metrics                            |                       |
| `test_prop_coefficient_dense_sparse_equivalence`                   | `n` identical dense rows equal one sparse site with multiplicity `n`               |                       |
| `test_prop_coefficient_additivity`                                 | Multi-site metrics equal sum of per-site metrics                                   |                       |
| `test_prop_coefficient_dense_finite_difference_derivative`         | Analytical first derivative matches central difference of `log_lh`                 |                       |
| `test_prop_coefficient_dense_finite_difference_second_derivative`  | Analytical Hessian matches second difference of `log_lh`                           | **1e-3** max_relative |
| `test_prop_coefficient_dense_hessian_matches_d1_finite_difference` | Analytical Hessian matches central difference of analytical first derivative       |                       |
| `test_hessian_stable_in_cancellation_regime`                       | Hessian preserves precision when posterior is concentrated on one eigenvalue class |                       |

### Support

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_dense/test_coefficient_extraction_dense_support.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_dense/test_coefficient_extraction_dense_support.rs)

Helper file. Provides `make_dense_seq_dis()`. No tests.

---

## Coefficient Extraction - Dense Invariants

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_dense/test_coefficient_extraction_dense_invariants.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_dense/test_coefficient_extraction_dense_invariants.rs)

**Impl:**

- [`packages/treetime/src/commands/optimize/optimize_dense.rs`](../../packages/treetime/src/commands/optimize/optimize_dense.rs)
- [`packages/treetime/src/commands/optimize/optimize_dense_eval.rs`](../../packages/treetime/src/commands/optimize/optimize_dense_eval.rs)
- [`packages/treetime/src/commands/optimize/optimize_sparse.rs`](../../packages/treetime/src/commands/optimize/optimize_sparse.rs)
- [`packages/treetime/src/commands/optimize/optimize_sparse_eval.rs`](../../packages/treetime/src/commands/optimize/optimize_sparse_eval.rs)

Parameterized invariant tests using rstest.

| Test                                                  | Purpose                                                            |
| ----------------------------------------------------- | ------------------------------------------------------------------ |
| `test_coefficient_invariant_nonneg_site_lh_at_zero`   | Coefficient sum non-negative at t=0 for valid probability messages |
| `test_coefficient_invariant_multiplicity_linearity`   | Sparse contribution with multiplicity m equals m times single-site |
| `test_coefficient_invariant_dense_sparse_equivalence` | m identical dense rows equal one sparse site with multiplicity m   |
| `test_coefficient_invariant_additivity`               | Two independent sites' log-likelihood equals sum of individual     |

---

## Coefficient Extraction - Sparse

**Directory:** [`test_coefficient_extraction_sparse/`](../../packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_sparse/)

### Site Tests

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_sparse/test_coefficient_extraction_sparse_site.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_sparse/test_coefficient_extraction_sparse_site.rs)

**Impl:** [`packages/treetime/src/commands/optimize/optimize_sparse.rs`](../../packages/treetime/src/commands/optimize/optimize_sparse.rs)

| Test                                                          | Purpose                         |
| ------------------------------------------------------------- | ------------------------------- |
| `test_site_contribution_stores_multiplicity_and_coefficients` | Store multiplicity and coeffs   |
| `test_site_contribution_unit_multiplicity`                    | Unit multiplicity for variables |

### Coefficients Tests

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_sparse/test_coefficient_extraction_sparse_coefficients.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_sparse/test_coefficient_extraction_sparse_coefficients.rs)

**Impl:**

- [`packages/treetime/src/commands/optimize/optimize_sparse.rs`](../../packages/treetime/src/commands/optimize/optimize_sparse.rs)
- [`packages/treetime/src/commands/optimize/optimize_sparse_eval.rs`](../../packages/treetime/src/commands/optimize/optimize_sparse_eval.rs)

| Test                                                       | Purpose                      |
| ---------------------------------------------------------- | ---------------------------- |
| `test_coefficients_computed_via_eigenvector_decomposition` | Eigenvector-based formula    |
| `test_matching_states_high_lh_at_zero`                     | High LH for matching states  |
| `test_mismatched_states_low_lh_at_zero`                    | Low LH for mismatched states |

### Derivatives Tests

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_sparse/test_coefficient_extraction_sparse_derivatives.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_sparse/test_coefficient_extraction_sparse_derivatives.rs)

**Impl:**

- [`packages/treetime/src/commands/optimize/optimize_dense.rs`](../../packages/treetime/src/commands/optimize/optimize_dense.rs)
- [`packages/treetime/src/commands/optimize/optimize_dense_eval.rs`](../../packages/treetime/src/commands/optimize/optimize_dense_eval.rs)
- [`packages/treetime/src/commands/optimize/optimize_sparse.rs`](../../packages/treetime/src/commands/optimize/optimize_sparse.rs)
- [`packages/treetime/src/commands/optimize/optimize_sparse_eval.rs`](../../packages/treetime/src/commands/optimize/optimize_sparse_eval.rs)

| Test                                                                 | Purpose                                         |
| -------------------------------------------------------------------- | ----------------------------------------------- |
| `test_evaluate_sparse_first_derivative_matches_numerical`            | Analytical vs numerical first derivative        |
| `test_evaluate_sparse_second_derivative_matches_numerical`           | Analytical vs numerical second derivative       |
| `test_evaluate_sparse_second_derivative_numerical_high_multiplicity` | Second derivative numerical at multiplicity=500 |
| `test_evaluate_sparse_derivatives_scale_with_multiplicity`           | All metrics scale linearly with multiplicity    |
| `test_evaluate_sparse_matches_dense_multiplicity_1`                  | Sparse (m=1) matches dense for all metrics      |

### Eigenvalues Tests

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_sparse/test_coefficient_extraction_sparse_eigenvalues.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_sparse/test_coefficient_extraction_sparse_eigenvalues.rs)

**Impl:**

- [`packages/treetime/src/commands/optimize/optimize_sparse.rs`](../../packages/treetime/src/commands/optimize/optimize_sparse.rs)
- [`packages/treetime/src/commands/optimize/optimize_sparse_eval.rs`](../../packages/treetime/src/commands/optimize/optimize_sparse_eval.rs)

| Test                                               | Purpose                          |
| -------------------------------------------------- | -------------------------------- |
| `test_eigenvalues_affect_branch_length_evaluation` | LH differs at different branches |
| `test_jc69_eigenvalues_structure`                  | One zero, three negative eigvals |

### Evaluate Tests

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_sparse/test_coefficient_extraction_sparse_evaluate.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_sparse/test_coefficient_extraction_sparse_evaluate.rs)

**Impl:**

- [`packages/treetime/src/commands/optimize/optimize_sparse.rs`](../../packages/treetime/src/commands/optimize/optimize_sparse.rs)
- [`packages/treetime/src/commands/optimize/optimize_sparse_eval.rs`](../../packages/treetime/src/commands/optimize/optimize_sparse_eval.rs)

| Test                                               | Purpose                     |
| -------------------------------------------------- | --------------------------- |
| `test_evaluate_sparse_single_site_multiplicity_1`  | Single site LH at zero      |
| `test_evaluate_sparse_single_site_multiplicity_10` | Multiplicity 10 LH at zero  |
| `test_evaluate_sparse_multiplicity_scales_log_lh`  | LH scales with multiplicity |
| `test_evaluate_sparse_multiple_sites_sum_log_lh`   | Multiple sites sum log-LH   |

### Mixed Sites Tests

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_sparse/test_coefficient_extraction_sparse_mixed_sites.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_sparse/test_coefficient_extraction_sparse_mixed_sites.rs)

**Impl:**

- [`packages/treetime/src/commands/optimize/optimize_sparse.rs`](../../packages/treetime/src/commands/optimize/optimize_sparse.rs)
- [`packages/treetime/src/commands/optimize/optimize_sparse_eval.rs`](../../packages/treetime/src/commands/optimize/optimize_sparse_eval.rs)

| Test                                               | Purpose                        |
| -------------------------------------------------- | ------------------------------ |
| `test_mixed_variable_and_fixed_sites`              | Variable + fixed sites compose |
| `test_fixed_sites_dominate_with_high_multiplicity` | High multiplicity dominates    |

### Partition Tests

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_sparse/test_coefficient_extraction_sparse_partition.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_coefficient_extraction_sparse/test_coefficient_extraction_sparse_partition.rs)

**Impl:** [`packages/treetime/src/commands/optimize/optimize_sparse.rs`](../../packages/treetime/src/commands/optimize/optimize_sparse.rs)

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

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_newton_convergence/test_newton_convergence_iteration.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_newton_convergence/test_newton_convergence_iteration.rs)

**Impl:** [`packages/treetime/src/commands/optimize/optimize_unified.rs`](../../packages/treetime/src/commands/optimize/optimize_unified.rs)

| Test                                            | Purpose                   |
| ----------------------------------------------- | ------------------------- |
| `test_newton_iteration_converges_within_bounds` | Convergence within bounds |
| `test_newton_iteration_respects_max_iter`       | Stops at max iterations   |

### Tolerance Tests

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_newton_convergence/test_newton_convergence_tolerance.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_newton_convergence/test_newton_convergence_tolerance.rs)

**Impl:**

- [`packages/treetime/src/commands/optimize/optimize_unified.rs`](../../packages/treetime/src/commands/optimize/optimize_unified.rs)
- [`packages/treetime/src/commands/optimize/method_newton.rs`](../../packages/treetime/src/commands/optimize/method_newton.rs)

| Test                                                     | Purpose                                    |
| -------------------------------------------------------- | ------------------------------------------ |
| `test_newton_tolerance_positive_branch_uses_relative`    | Relative tolerance dominates for large BL  |
| `test_newton_tolerance_zero_branch_uses_absolute_floor`  | Absolute floor prevents degenerate zero    |
| `test_newton_tolerance_small_branch_uses_absolute_floor` | Absolute floor for very small BL           |
| `test_newton_tolerance_crossover`                        | Relative = absolute at crossover point     |
| `test_newton_tolerance_always_positive`                  | Positive and finite for all BL values      |
| `test_newton_zero_start_converges_efficiently`           | Zero-start does not exhaust all iterations |

### Evaluate Mixed Tests

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_newton_convergence/test_newton_convergence_evaluate_mixed.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_newton_convergence/test_newton_convergence_evaluate_mixed.rs)

**Impl:** [`packages/treetime/src/commands/optimize/optimize_unified.rs`](../../packages/treetime/src/commands/optimize/optimize_unified.rs)

| Test                                          | Purpose                       |
| --------------------------------------------- | ----------------------------- |
| `test_evaluate_mixed_returns_finite_values`   | Finite values for valid input |
| `test_evaluate_mixed_log_lh_negative`         | Log-LH is non-positive        |
| `test_evaluate_mixed_multiple_branch_lengths` | Finite across branch lengths  |

### Support

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_newton_convergence/test_newton_convergence_support.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_newton_convergence/test_newton_convergence_support.rs)

Helper file. Provides `make_dense_contribution()`. No tests.

---

## Grid Search

**Directory:** [`test_grid_search/`](../../packages/treetime/src/commands/optimize/__tests__/test_grid_search/)

### Basic Tests

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_grid_search/test_grid_search_basic.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_grid_search/test_grid_search_basic.rs)

**Impl:** [`packages/treetime/src/commands/optimize/optimize_unified.rs`](../../packages/treetime/src/commands/optimize/optimize_unified.rs)

| Test                                        | Purpose                     |
| ------------------------------------------- | --------------------------- |
| `test_grid_search_finds_maximum_log_lh`     | Finds likelihood maximum    |
| `test_grid_search_result_in_range`          | Result within search bounds |
| `test_grid_search_with_small_branch_length` | Works with small branches   |
| `test_grid_search_with_large_branch_length` | Works with large branches   |

### Edge Cases Tests

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_grid_search/test_grid_search_edge_cases.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_grid_search/test_grid_search_edge_cases.rs)

**Impl:** [`packages/treetime/src/commands/optimize/optimize_unified.rs`](../../packages/treetime/src/commands/optimize/optimize_unified.rs)

| Test                                                | Purpose                      |
| --------------------------------------------------- | ---------------------------- |
| `test_grid_search_one_mutation_dominates_range`     | Large one_mutation vs branch |
| `test_grid_search_preserves_positive_branch_length` | Always returns positive      |

### Partitions Tests

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_grid_search/test_grid_search_partitions.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_grid_search/test_grid_search_partitions.rs)

**Impl:** [`packages/treetime/src/commands/optimize/optimize_unified.rs`](../../packages/treetime/src/commands/optimize/optimize_unified.rs)

| Test                                    | Purpose                         |
| --------------------------------------- | ------------------------------- |
| `test_grid_search_multiple_partitions`  | Combines multiple contributions |
| `test_grid_search_uniform_coefficients` | Handles uniform coefficients    |

### Coverage Tests

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_grid_search/test_grid_search_coverage.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_grid_search/test_grid_search_coverage.rs)

**Impl:** [`packages/treetime/src/commands/optimize/optimize_unified.rs`](../../packages/treetime/src/commands/optimize/optimize_unified.rs)

| Test                                                                  | Purpose                                                  |
| --------------------------------------------------------------------- | -------------------------------------------------------- |
| `test_grid_search_branch_lengths_zero_branch_covers_full_range`       | Grid extends to 0.5 subs/site when branch_length=0       |
| `test_grid_search_branch_lengths_large_branch_extends_beyond_minimum` | Proportional bound dominates for large branches          |
| `test_grid_search_branch_lengths_is_log_spaced`                       | Grid is log-spaced with constant ratio                   |
| `test_grid_search_branch_lengths_has_correct_point_count`             | Grid always has 100 points                               |
| `test_grid_search_reaches_extended_range_from_zero`                   | Monotonic likelihood finds optimum near upper bound      |
| `test_grid_search_branch_lengths_coverage_invariant`                  | Grid covers [one_mutation, 0.5] for all parameter combos |

### Support

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_grid_search/test_grid_search_support.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_grid_search/test_grid_search_support.rs)

Helper file. Provides `make_dense_contribution()` and `grid_search()`. No tests.

---

## Dense/Sparse Equivalence

> **Note**: All tests in this directory are parameterized across 6 `BranchOptMethod` variants: Newton, NewtonSqrt, NewtonLog, Brent, BrentSqrt, BrentLog.

**Directory:** [`test_dense_sparse_equivalence/`](../../packages/treetime/src/commands/optimize/__tests__/test_dense_sparse_equivalence/)

### Initial Tests

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_initial.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_initial.rs)

**Impl:** [`packages/treetime/src/commands/optimize/optimize_unified.rs`](../../packages/treetime/src/commands/optimize/optimize_unified.rs)

| Test                                                          | Purpose                      |
| ------------------------------------------------------------- | ---------------------------- |
| `test_dense_sparse_initial_log_lh_equivalence`                | Initial LH matches           |
| `test_dense_sparse_initial_log_lh_equivalence_with_mutations` | Initial LH matches with muts |

### Convergence Tests

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_convergence.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_convergence.rs)

**Impl:** [`packages/treetime/src/commands/optimize/optimize_unified.rs`](../../packages/treetime/src/commands/optimize/optimize_unified.rs)

| Test                                 | Purpose               |
| ------------------------------------ | --------------------- |
| `test_dense_optimization_converges`  | Dense mode converges  |
| `test_sparse_optimization_converges` | Sparse mode converges |

### Bounds Tests

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_bounds.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_bounds.rs)

**Impl:** [`packages/treetime/src/commands/optimize/optimize_unified.rs`](../../packages/treetime/src/commands/optimize/optimize_unified.rs)

| Test                                                             | Purpose                    |
| ---------------------------------------------------------------- | -------------------------- |
| `test_dense_sparse_log_lh_bounded_difference_after_optimization` | LH difference is bounded   |
| `test_dense_sparse_branch_lengths_bounded_difference`            | Branch lengths are similar |

### Validity Tests

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_validity.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_validity.rs)

**Impl:** [`packages/treetime/src/commands/optimize/optimize_unified.rs`](../../packages/treetime/src/commands/optimize/optimize_unified.rs)

| Test                                              | Purpose                       |
| ------------------------------------------------- | ----------------------------- |
| `test_dense_optimization_produces_valid_results`  | Dense produces valid results  |
| `test_sparse_optimization_produces_valid_results` | Sparse produces valid results |

### Support

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_support.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_dense_sparse_equivalence/test_dense_sparse_equivalence_support.rs)

Helper file. Provides `setup_dense_only()`, `setup_sparse_only()`, `get_branch_lengths()`, `gap_free_alignment()`. No tests.

---

## Convergence Control

> **Note**: All tests in this directory are parameterized across 6 `BranchOptMethod` variants: Newton, NewtonSqrt, NewtonLog, Brent, BrentSqrt, BrentLog.

**Directory:** [`test_convergence/`](../../packages/treetime/src/commands/optimize/__tests__/test_convergence/)

### Iterations Tests

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_convergence/test_convergence_iterations.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_convergence/test_convergence_iterations.rs)

**Impl:** [`packages/treetime/src/commands/optimize/optimize_unified.rs`](../../packages/treetime/src/commands/optimize/optimize_unified.rs)

| Test                                                 | Purpose                                                   | Notes                         |
| ---------------------------------------------------- | --------------------------------------------------------- | ----------------------------- |
| `test_optimization_converges_within_iterations`      | Bounded oscillation: undamped LH stays within tight range |                               |
| `test_optimization_improves_or_maintains_likelihood` | Strict non-regression: final LH >= initial LH             |                               |
| `test_optimization_produces_valid_branch_lengths`    | Branch lengths valid after opt                            | **1e-2** cross-method LH diff |

### Edge Cases Tests

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_convergence/test_convergence_edge_cases.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_convergence/test_convergence_edge_cases.rs)

**Impl:** [`packages/treetime/src/commands/optimize/optimize_unified.rs`](../../packages/treetime/src/commands/optimize/optimize_unified.rs)

| Test                                            | Purpose                      |
| ----------------------------------------------- | ---------------------------- |
| `test_optimization_handles_zero_branch_lengths` | Handles zero-length branches |
| `test_optimization_handles_very_short_branches` | Handles very short branches  |
| `test_optimization_handles_long_branches`       | Handles long branches        |

### Idempotence Tests

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_convergence/test_convergence_idempotence.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_convergence/test_convergence_idempotence.rs)

**Impl:** [`packages/treetime/src/commands/optimize/optimize_unified.rs`](../../packages/treetime/src/commands/optimize/optimize_unified.rs)

| Test                                                    | Purpose                        |
| ------------------------------------------------------- | ------------------------------ |
| `test_optimization_converges_with_valid_branch_lengths` | Branch lengths valid each iter |
| `test_second_optimization_produces_same_likelihood`     | Deterministic across runs      |

### Support

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_convergence/test_convergence_support.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_convergence/test_convergence_support.rs)

Helper file. Provides `simple_alignment()`, `setup_partitions()`, `compute_total_lh()`. No tests.

---

## Outer-Loop Damping

> **Note**: `test_damped_optimization_converges` and `test_damped_optimization_does_not_regress` are parameterized across 6 `BranchOptMethod` variants.

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_damping.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_damping.rs)

**Impl:** [`packages/treetime/src/commands/optimize/run.rs`](../../packages/treetime/src/commands/optimize/run.rs)

| Test                                                     | Purpose                                                                           |
| -------------------------------------------------------- | --------------------------------------------------------------------------------- |
| `test_save_branch_lengths_captures_all_edges`            | Verifies save_branch_lengths reads all edges                                      |
| `test_apply_damping_zero_is_noop`                        | damping=0.0 leaves optimized values unchanged                                     |
| `test_apply_damping_weights_match_v0` (5 cases)          | Blend weights match v0 formula at iterations 0-9                                  |
| `test_apply_damping_blends_correctly`                    | Verifies bl = new*w_new + old*w_old                                               |
| `test_apply_damping_new_weight_increases_with_iteration` | New weight monotonically increases over iterations                                |
| `test_damped_optimization_converges`                     | Production convergence criterion, sign-flip oscillation detection, tail stability |
| `test_damped_optimization_does_not_regress`              | Strict non-regression: final LH >= initial LH                                     |

---

## v0 Parity and Damped vs Undamped (Golden Master)

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_gm_optimize.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_gm_optimize.rs)

**Impl:**

- [`packages/treetime/src/commands/optimize/optimize_unified.rs`](../../packages/treetime/src/commands/optimize/optimize_unified.rs)
- [`packages/treetime/src/commands/optimize/run.rs`](../../packages/treetime/src/commands/optimize/run.rs)

| Test                                  | Purpose                                                        | Notes                                          |
| ------------------------------------- | -------------------------------------------------------------- | ---------------------------------------------- |
| `test_gm_optimize`                    | v0/v1 total branch length parity within 5% (flu/h3n2/20, JC69) | **ignored**: per-branch divergence exceeds 10% |
| `test_gm_optimize_damped_vs_undamped` | Damping does not increase oscillation on real dataset          |                                                |

Fixtures in `__fixtures__/`:

- [`gm_optimize_inputs.json`](../../packages/treetime/src/commands/optimize/__tests__/__fixtures__/gm_optimize_inputs.json) - shared test parameters
- [`gm_optimize_outputs.json`](../../packages/treetime/src/commands/optimize/__tests__/__fixtures__/gm_optimize_outputs.json) - v0 reference captured by `gm_optimize_capture`
- `gm_optimize_capture` - Python capture script for v0 TreeAnc.optimize_tree_marginal

---

## Indel Objective and Zero-BL Handling

**File:** [`test_optimize_indel.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_optimize_indel.rs)

Focused regression coverage for the Poisson indel term, zero-branch-length bootstrapping, and indel-aware optimization behavior.

| Test                                                     | Purpose                                                                  |
| -------------------------------------------------------- | ------------------------------------------------------------------------ |
| `test_optimize_indel_estimate_rate_no_indels`            | Global indel rate is zero when no edges carry indels                     |
| `test_optimize_indel_estimate_rate_with_indels`          | Global indel rate equals total indels divided by total branch length     |
| `test_optimize_indel_total_log_lh_matches_manual_sum`    | Tree-level Poisson indel objective matches manual per-edge summation     |
| `test_optimize_indel_initial_guess_nonzero_with_indels`  | Indel-only signal produces positive initial branch length                |
| `test_optimize_indel_initial_guess_zero_bl_tree_with_indels` | Zero-BL trees bootstrap away from zero when indels are present       |
| `test_optimize_indel_run_optimize_nonzero_with_indels`   | All 6 optimize methods produce positive finite BL on indel-bearing edge  |
| `test_optimize_indel_zero_bl_pipeline_escapes_zero`      | Initial guess plus optimize escapes zero in the sparse production path   |
| `test_optimize_indel_grid_zero_comparison_rejects_zero_with_indels` | Zero-boundary shortcut rejects zero when Poisson term is active |
| `test_optimize_indel_grid_zero_comparison_allows_zero_without_indels` | Zero-boundary shortcut preserves substitution-only behavior      |

---

## Run Optimize Loop Contract

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_run_optimize_loop.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_run_optimize_loop.rs)

**Impl:** [`packages/treetime/src/commands/optimize/run.rs`](../../packages/treetime/src/commands/optimize/run.rs)

| Test                                              | Purpose                                                  |
| ------------------------------------------------- | -------------------------------------------------------- |
| `test_run_optimize_loop_records_lh_history`                         | `lh_history` has one entry per executed iteration                    |
| `test_run_optimize_loop_records_joint_likelihood_with_sparse_indels` | `lh_history` records substitution + indel objective on indel runs |
| `test_run_optimize_loop_breaks_on_convergence`                      | Breaks and records `stopped_at` with `Converged`                     |
| `test_run_optimize_loop_zero_max_iter_is_noop`                      | `max_iter = 0` runs the body zero times                              |
| `test_run_optimize_loop_all_likelihoods_finite`                     | All likelihoods finite (guards against forward-pass NaN)             |
| `test_run_optimize_loop_improves_likelihood`                        | Damped loop improves likelihood from initial state                   |

---

## Convergence Conditions

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_convergence_conditions.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_convergence_conditions.rs)

**Impl:** [`packages/treetime/src/commands/optimize/run.rs`](../../packages/treetime/src/commands/optimize/run.rs)

| Test                                                                    | Purpose                                                |
| ----------------------------------------------------------------------- | ------------------------------------------------------ |
| `test_convergence_conditions_damping_floor_at_high_iteration` (3 cases) | DAMPING_FLOOR weight at iterations 100/500/1000        |
| `test_convergence_conditions_damping_uses_exponential_below_crossover`  | Exponential decay used below crossover (iteration 5)   |
| `test_convergence_conditions_restore_branch_lengths_roundtrip`          | Save/modify/restore/verify identical                   |
| `test_convergence_conditions_converged_reason`                          | Converged or Oscillating on damped toy tree            |
| `test_convergence_conditions_worsened_reverts_to_best`                  | Worsened fires on undamped toy tree, trigger LH < best |
| `test_convergence_conditions_oscillation_detection`                     | Early stop with large dp                               |
| `test_convergence_conditions_exhausts_max_iter`                         | max_iter=2 exhausted without stopping                  |
| `test_convergence_conditions_dense_only_converges`                      | Dense-only regression check                            |

---

## Convergence on Real Datasets

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_convergence_sc2.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_convergence_sc2.rs)

**Impl:**

- [`packages/treetime/src/commands/optimize/optimize_unified.rs`](../../packages/treetime/src/commands/optimize/optimize_unified.rs)
- [`packages/treetime/src/commands/optimize/run.rs`](../../packages/treetime/src/commands/optimize/run.rs)

| Test                                                | Purpose                                              |
| --------------------------------------------------- | ---------------------------------------------------- |
| `test_convergence_sc2_sparse_converges_on_sc2_2844` | Regression: sc2/2844 sparse convergence (slow, ~60s) |
| `test_convergence_sc2_flu_h3n2_20_converges`        | flu/h3n2/20 sparse convergence within 10 iterations  |

---

## Optimization Metrics

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_optimization_metrics.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_optimization_metrics.rs)

**Impl:** [`packages/treetime/src/commands/optimize/optimize_unified.rs`](../../packages/treetime/src/commands/optimize/optimize_unified.rs)

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

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_is_zero_branch_optimal.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_is_zero_branch_optimal.rs)

**Impl:**

- [`packages/treetime/src/commands/optimize/optimize_unified.rs`](../../packages/treetime/src/commands/optimize/optimize_unified.rs)
- [`packages/treetime/src/commands/optimize/optimize_dense.rs`](../../packages/treetime/src/commands/optimize/optimize_dense.rs)
- [`packages/treetime/src/commands/optimize/optimize_sparse.rs`](../../packages/treetime/src/commands/optimize/optimize_sparse.rs)

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

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_initial_guess_mode.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_initial_guess_mode.rs)

**Impl:**

- [`packages/treetime/src/commands/optimize/optimize_unified.rs`](../../packages/treetime/src/commands/optimize/optimize_unified.rs)
- [`packages/treetime/src/commands/optimize/run.rs`](../../packages/treetime/src/commands/optimize/run.rs)
- [`packages/treetime/src/commands/optimize/args.rs`](../../packages/treetime/src/commands/optimize/args.rs)

| Test                                                            | Purpose                                                                 |
| --------------------------------------------------------------- | ----------------------------------------------------------------------- |
| `test_initial_guess_mode_default_is_auto`                       | Default variant is Auto                                                 |
| `test_initial_guess_mode_detects_nan_from_newick`               | Detects NaN branch lengths from bio crate newick parser                 |
| `test_initial_guess_mode_detects_explicit_nan`                  | Detects explicit `Some(NaN)` as missing                                 |
| `test_initial_guess_mode_no_missing_when_all_finite`            | Reports no missing when all edges have finite lengths                   |
| `test_initial_guess_mode_auto_preserves_valid_lengths`          | Auto mode does not overwrite valid finite branch lengths                |
| `test_initial_guess_mode_auto_fills_nan_from_newick`            | Auto mode fills NaN edges with finite non-negative values               |
| `test_initial_guess_mode_auto_fills_only_missing_edges`         | Auto mode fills one NaN edge, preserves all others exactly              |
| `test_initial_guess_mode_always_overwrites_all`                 | Always mode overwrites existing branch lengths                          |
| `test_initial_guess_mode_never_accepts_complete_tree`           | Never mode accepts tree with all finite branch lengths                  |
| `test_initial_guess_mode_never_rejects_nan_tree`                | Never mode detects NaN tree as having missing branch lengths            |
| `test_initial_guess_mode_never_accepts_zero_bl_without_indels`  | Never mode accepts all-zero branch lengths when no indels are present   |
| `test_initial_guess_mode_never_rejects_zero_bl_with_indels`     | Never mode rejects zero branch length on an indel-bearing edge          |
| `test_initial_guess_mode_never_accepts_positive_bl_with_indels` | Never mode accepts positive branch lengths on indel-bearing edges       |
| `test_any_indel_edge_has_zero_bl_false_without_indels`          | Helper returns false when no indels are present                         |
| `test_any_indel_edge_has_zero_bl_false_with_positive_bl`        | Helper returns false when indel-bearing edge has positive branch length |
| `test_any_indel_edge_has_zero_bl_true_with_indel_and_zero_bl`   | Helper returns true when an indel-bearing edge has zero branch length   |

---

## Initial Guess Indel Zero-BL

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_initial_guess_indel_zero_bl.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_initial_guess_indel_zero_bl.rs)

**Impl:** [`packages/treetime/src/commands/optimize/optimize_unified.rs`](../../packages/treetime/src/commands/optimize/optimize_unified.rs)

| Test                                                       | Purpose                                            |
| ---------------------------------------------------------- | -------------------------------------------------- |
| `test_initial_guess_auto_preserves_zero_bl_without_indels` | Zero BL preserved when no indels present           |
| `test_initial_guess_auto_overrides_zero_bl_with_indels`    | Zero BL overridden to positive when indels present |

---

## Initial Guess GTR Messages

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_initial_guess_gtr_messages.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_initial_guess_gtr_messages.rs)

**Impl:** [`packages/treetime/src/commands/optimize/optimize_unified.rs`](../../packages/treetime/src/commands/optimize/optimize_unified.rs)

| Test                                             | Purpose                                                       |
| ------------------------------------------------ | ------------------------------------------------------------- |
| `test_stale_jc69_messages_bias_initial_guess`    | Stale JC69 messages produce different branch lengths than F81 |
| `test_initial_guess_idempotent_after_gtr_update` | Identical initialization sequences produce identical results  |

---

## Initial Guess Formula

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_initial_guess_formula.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_initial_guess_formula.rs)

**Impl:**

- [`packages/treetime/src/commands/optimize/optimize_unified.rs`](../../packages/treetime/src/commands/optimize/optimize_unified.rs)
- [`packages/treetime/src/commands/optimize/partition_ops.rs`](../../packages/treetime/src/commands/optimize/partition_ops.rs)

| Test                                                                                 | Purpose                                                                                  |
| ------------------------------------------------------------------------------------ | ---------------------------------------------------------------------------------------- |
| `test_initial_guess_formula_sparse`                                                  | Sparse branch length equals substitution count over effective length                     |
| `test_initial_guess_formula_dense`                                                   | Dense branch length equals substitution count over effective length                      |
| `test_initial_guess_dense_sparse_ambiguous_r_reference_state_consistency`            | Dense and sparse initial branch lengths match on partial ambiguity                       |
| `test_optimize_contribution_dense_sparse_ambiguous_r_value_and_gradient_consistency` | Dense and sparse optimize contributions agree in value and gradient on partial ambiguity |

---

## Initial Guess Gaps

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_initial_guess_gaps.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_initial_guess_gaps.rs)

**Impl:** [`packages/treetime/src/commands/optimize/optimize_unified.rs`](../../packages/treetime/src/commands/optimize/optimize_unified.rs)

| Test                                           | Purpose                                                                     |
| ---------------------------------------------- | --------------------------------------------------------------------------- |
| `test_sparse_effective_length_no_gaps`         | Sparse effective length equals full length with no gaps                     |
| `test_dense_effective_length_no_gaps`          | Dense effective length equals full length with no gaps                      |
| `test_sparse_effective_length_shared_gaps`     | Shared gaps reduce sparse effective length by gap count                     |
| `test_dense_effective_length_shared_gaps`      | Shared gaps reduce dense effective length by gap count                      |
| `test_sparse_effective_length_one_leaf_gapped` | One-leaf gaps reduce effective length on that edge only                     |
| `test_dense_edge_subs_excludes_gap_positions`  | Dense edge_subs excludes gap positions from substitutions                   |
| `test_initial_guess_sparse_gap_adjusted_rate`  | Initial guess adjusts branch length rate proportionally to effective length |

---

## Dense Edge Subs

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_dense_edge_subs.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_dense_edge_subs.rs)

**Impl:** [`packages/treetime/src/representation/partition/marginal_dense.rs`](../../packages/treetime/src/representation/partition/marginal_dense.rs)

| Test                                                                 | Purpose                                                                         |
| -------------------------------------------------------------------- | ------------------------------------------------------------------------------- |
| `test_dense_edge_subs_no_false_mutation_from_uniform_outgroup`       | Uniform outgroup message does not create false substitutions                    |
| `test_dense_edge_subs_detects_real_mutation_hidden_by_edge_messages` | Edge messages disagreeing with posteriors do not mask real substitutions        |
| `test_dense_edge_subs_match_reconstructed_branch_differences`        | Dense edge_subs match parent-child MAP sequence differences after full marginal |
| `test_dense_edge_subs_excludes_gap_positions_with_posteriors`        | Gap positions excluded from dense edge_subs even when posteriors differ         |
| `test_dense_edge_subs_is_canonical_filter_present`                   | is_canonical filter validates only canonical states appear in substitutions     |

---

## Eval Zero-Branch Mismatch

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_eval_zero_branch_mismatch.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_eval_zero_branch_mismatch.rs)

**Impl:**

- [`packages/treetime/src/commands/optimize/optimize_unified.rs`](../../packages/treetime/src/commands/optimize/optimize_unified.rs)
- [`packages/treetime/src/commands/optimize/run.rs`](../../packages/treetime/src/commands/optimize/run.rs)

| Test                                    | Purpose                                                                                 |
| --------------------------------------- | --------------------------------------------------------------------------------------- |
| `test_eval_zero_branch_mismatch_no_nan` | run_optimize_mixed does not produce -inf/NaN with zero BL and mismatched certain states |

---

## Topology Cleanup in Loop

> **Note**: `test_optimize_loop_with_topology_cleanup_sparse`, `test_optimize_loop_no_collapse_when_branches_nonzero`, and `test_optimize_loop_with_topology_cleanup_dense` are parameterized across 6 `BranchOptMethod` variants.

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_topology_cleanup.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_topology_cleanup.rs)

**Impl:**

- [`packages/treetime/src/commands/optimize/optimize_unified.rs`](../../packages/treetime/src/commands/optimize/optimize_unified.rs)
- [`packages/treetime/src/commands/optimize/run.rs`](../../packages/treetime/src/commands/optimize/run.rs)

| Test                                                               | Purpose                                                    |
| ------------------------------------------------------------------ | ---------------------------------------------------------- |
| `test_optimize_find_zero_optimal_internal_edges_empty_graph`       | Empty graph returns no edges                               |
| `test_optimize_find_zero_optimal_internal_edges_no_zero_edges`     | Non-zero branches return no edges                          |
| `test_optimize_find_zero_optimal_internal_edges_skips_leaves`      | Leaf edges with bl=0 are not collected                     |
| `test_optimize_find_zero_optimal_internal_edges_collects_internal` | Internal edges with bl=0 are collected                     |
| `test_optimize_find_zero_optimal_internal_edges_multiple`          | Multiple zero-optimal internal edges collected             |
| `test_optimize_prune_and_merge_empty_list`                         | Empty zero-edge list is a noop                             |
| `test_optimize_prune_and_merge_collapses_and_merges`               | Collapse + shared mutation merge in one pass               |
| `test_optimize_loop_with_topology_cleanup_sparse`                  | Full loop collapses zero-optimal branches (identical seqs) |
| `test_optimize_loop_no_collapse_when_branches_nonzero`             | No collapse when all branches carry signal                 |
| `test_optimize_merge_then_marginal_finite_likelihood`              | Merge + update_marginal produces finite log-likelihood     |
| `test_optimize_cascading_collapse_parent_child_both_zero`          | Parent-child both zero-optimal, guard handles removal      |
| `test_optimize_loop_with_topology_cleanup_dense`                   | Dense-mode full loop with topology cleanup                 |

Direct coverage for `collapse_edge()` lives in [Representation Tests: Topology Cleanup / Edge Collapse](representation.md#topology-cleanup--edge-collapse) (8 tests).

---

## Indel Contribution (Inline)

**Test:** [`packages/treetime/src/commands/optimize/optimize_indel.rs`](../../packages/treetime/src/commands/optimize/optimize_indel.rs) (inline `#[cfg(test)]`)

**Impl:** [`packages/treetime/src/commands/optimize/optimize_indel.rs`](../../packages/treetime/src/commands/optimize/optimize_indel.rs)

| Test                                                        | Purpose                                                           |
| ----------------------------------------------------------- | ----------------------------------------------------------------- |
| `test_optimize_indel_poisson_zero_rate`                     | Poisson metrics are all zero when mu=0                            |
| `test_optimize_indel_poisson_zero_indels`                   | k=0 gives log_lh=-mu\*t, derivative=-mu, second_derivative=0      |
| `test_optimize_indel_poisson_log_lh_value`                  | Specific log-likelihood value for k=2, mu=10, t=0.1               |
| `test_optimize_indel_poisson_derivative`                    | Analytical derivative and second derivative match expected values |
| `test_optimize_indel_poisson_mle_at_optimum`                | Derivative is zero at MLE t=k/mu                                  |
| `test_optimize_indel_poisson_derivative_positive_near_zero` | Derivative is large positive near t=0 for k>0                     |
| `test_optimize_indel_poisson_second_derivative_negative`    | Second derivative negative for k>0 (log-concave)                  |
| `test_optimize_indel_statrs_ln_factorial`                   | statrs ln_factorial agrees with direct computation                |

---

## Indel Contribution (Integration)

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_optimize_indel.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_optimize_indel.rs)

**Impl:**

- [`packages/treetime/src/commands/optimize/optimize_indel.rs`](../../packages/treetime/src/commands/optimize/optimize_indel.rs)
- [`packages/treetime/src/commands/optimize/optimize_unified.rs`](../../packages/treetime/src/commands/optimize/optimize_unified.rs)
- [`packages/treetime/src/commands/optimize/optimize_dense.rs`](../../packages/treetime/src/commands/optimize/optimize_dense.rs)
- [`packages/treetime/src/commands/optimize/run.rs`](../../packages/treetime/src/commands/optimize/run.rs)

| Test                                                                  | Purpose                                                                        |
| --------------------------------------------------------------------- | ------------------------------------------------------------------------------ |
| `test_optimize_indel_estimate_rate_no_indels`                         | estimate_indel_rate returns 0 with no indels                                   |
| `test_optimize_indel_estimate_rate_with_indels`                       | estimate_indel_rate returns total_indels / total_branch_length                 |
| `test_optimize_indel_initial_guess_nonzero_with_indels`               | initial_guess assigns positive BL when indels present on identical sequences   |
| `test_optimize_indel_initial_guess_zero_bl_tree_with_indels`          | initial_guess bootstraps positive BL on zero-BL tree with indels               |
| `test_optimize_indel_run_optimize_nonzero_with_indels` (6 methods)    | run_optimize_mixed assigns positive BL with indels for each optimizer          |
| `test_optimize_indel_zero_bl_pipeline_escapes_zero` (6 methods)       | Full pipeline escapes zero BL with indels for each optimizer                   |
| `test_optimize_indel_poisson_concavity` (5 cases)                     | Poisson second derivative negative for k>0                                     |
| `test_optimize_indel_poisson_mle_derivative_zero` (5 cases)           | Poisson derivative zero at MLE t=k/mu                                          |
| `test_optimize_indel_poisson_mle_is_maximum` (6 cases)                | Poisson log-likelihood at MLE exceeds any other point                          |
| `test_optimize_indel_poisson_numerical_derivative`                    | Numerical derivative matches analytical derivative                             |
| `test_optimize_indel_zero_branch_no_indels_unchanged`                 | is_zero_branch_optimal returns true with no indels and negative sub derivative |
| `test_optimize_indel_grid_zero_comparison_rejects_zero_with_indels`   | Grid search rejects zero when indels present                                   |
| `test_optimize_indel_grid_zero_comparison_allows_zero_without_indels` | Grid search allows zero without indels when subs prefer it                     |
| `test_optimize_indel_newton_converges_to_poisson_mle` (4 cases)       | Newton step from nearby point converges toward Poisson MLE                     |
| `test_optimize_indel_min_branch_length_clamping` (6 methods)          | Optimizer clamps BL above zero when indels present                             |
| `test_optimize_indel_evaluate_with_indels_shifts_log_lh` (4 cases)    | Indel term shifts combined log-likelihood from substitution-only value         |
| `test_prop_optimize_indel_concavity`                                  | Property: second derivative always negative for k>0                            |
| `test_prop_optimize_indel_mle_derivative`                             | Property: derivative is zero at MLE                                            |
| `test_prop_optimize_indel_derivative_positive_near_zero`              | Property: derivative positive near t=0 for k>0                                 |
| `test_prop_optimize_indel_numerical_derivative`                       | Property: numerical derivative matches analytical                              |

---

## Optimization Method

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_optimize_method.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_optimize_method.rs)

**Impl:**

- [`packages/treetime/src/commands/optimize/method_newton.rs`](../../packages/treetime/src/commands/optimize/method_newton.rs)
- [`packages/treetime/src/commands/optimize/method_brent.rs`](../../packages/treetime/src/commands/optimize/method_brent.rs)
- [`packages/treetime/src/commands/optimize/optimize_unified.rs`](../../packages/treetime/src/commands/optimize/optimize_unified.rs)
- [`packages/treetime/src/commands/optimize/optimize_dense.rs`](../../packages/treetime/src/commands/optimize/optimize_dense.rs)
- [`packages/treetime/src/commands/optimize/optimize_indel.rs`](../../packages/treetime/src/commands/optimize/optimize_indel.rs)

Tests the 6 per-edge branch length optimization methods (Newton, NewtonSqrt, NewtonLog, Brent, BrentSqrt, BrentLog) against 5 verification criteria:

- **C1: Local optimality**: log-likelihood at optimum exceeds neighbors at 3 delta scales (0.1%, 1%, 10%)
- **C2: Stationarity** (Newton-specific): implied Newton step at optimum is below tolerance
- **C3: Cross-method log-likelihood agreement**: all 6 methods produce log-likelihood within 1e-3
- **C4: Bracket validity** (Brent-specific): optimum log-likelihood exceeds bracket endpoints
- **C5: Cross-conditioning ordering**: on indel-bearing edges: `lh_newton_log >= lh_newton_sqrt >= lh_newton`

### Chain Rule Tests (sqrt and log transforms)

| Test                                                                        | Purpose                                             | Notes            |
| --------------------------------------------------------------------------- | --------------------------------------------------- | ---------------- |
| `test_optimize_method_chain_rule_at_zero`                                   | Chain rule sqrt at s=0: ds=0, d2s=2\*dl_dt          |                  |
| `test_optimize_method_chain_rule_analytical`                                | Chain rule sqrt at known analytical values          |                  |
| `test_optimize_method_chain_rule_log_analytical`                            | Chain rule log at known analytical values           |                  |
| `test_optimize_method_chain_rule_log_small_t`                               | Chain rule log approaches zero for small t          |                  |
| `test_optimize_method_chain_rule_numerical_first_derivative` (4 cases)      | sqrt first derivative matches numerical difference  |                  |
| `test_optimize_method_chain_rule_numerical_second_derivative` (4 cases)     | sqrt second derivative matches numerical difference | **1e-2** epsilon |
| `test_optimize_method_chain_rule_log_numerical_first_derivative` (4 cases)  | log first derivative matches numerical difference   |                  |
| `test_optimize_method_chain_rule_log_numerical_second_derivative` (4 cases) | log second derivative matches numerical difference  | **1e-2** epsilon |

### Method Equivalence and Optimality (C1, C3)

| Test                                                                  | Purpose                                           | Notes    |
| --------------------------------------------------------------------- | ------------------------------------------------- | -------- |
| `test_optimize_method_equivalence_no_indels` (6 cases)                | All 6 methods produce finite non-negative BLs     |          |
| `test_optimize_method_local_optimality` (6 cases)                     | C1: LH at optimum exceeds neighbors (all methods) |          |
| `test_optimize_method_cross_method_lh_agreement` (3 cases)            | C3: NewtonSqrt and Brent LH agree within 1e-3     | **1e-3** |
| `test_optimize_method_cross_method_lh_agreement_newton_log` (3 cases) | C3: NewtonLog and Brent LH agree within 1e-3      | **1e-3** |
| `test_optimize_method_cross_method_all_six_lh_agreement` (3 cases)    | C3: All 6 methods agree within 1e-3 vs BrentSqrt  | **1e-3** |

### Newton Stationarity and Conditioning (C2, C5)

| Test                                                                | Purpose                                              |
| ------------------------------------------------------------------- | ---------------------------------------------------- |
| `test_optimize_method_stationarity` (3 cases)                       | C2: Implied Newton step below tolerance (3 variants) |
| `test_optimize_method_newton_sqrt_improves_over_newton`             | C5: NewtonSqrt LH >= Newton on indel case            |
| `test_optimize_method_newton_log_improves_over_newton`              | C5: NewtonLog LH >= Newton on indel case             |
| `test_optimize_method_newton_cross_conditioning_ordering` (3 cases) | C5: lh_log >= lh_sqrt >= lh_t (all indel counts)     |

### Brent Bracket and Transform Validity (C4)

| Test                                                             | Purpose                                             |
| ---------------------------------------------------------------- | --------------------------------------------------- |
| `test_optimize_method_brent_bracket_validity` (3 cases)          | C4: Optimum LH exceeds bracket endpoints (3 Brents) |
| `test_optimize_method_brent_cross_parameterization_lh_agreement` | Brent-t, Brent-sqrt, Brent-log agree within 1e-3    |
| `test_optimize_method_brent_sqrt_transform_round_trip`           | sqrt transform produces local optimum in t-space    |
| `test_optimize_method_brent_log_transform_round_trip`            | log transform produces local optimum in t-space     |

### Indel Robustness

| Test                                                              | Purpose                                            |
| ----------------------------------------------------------------- | -------------------------------------------------- |
| `test_optimize_method_brent_positive_with_indels` (9 cases)       | All 3 Brent variants x 3 indel counts: positive BL |
| `test_optimize_method_newton_positive_with_indels` (3 cases)      | Newton-t produces positive finite BL with indels   |
| `test_optimize_method_newton_sqrt_positive_with_indels` (3 cases) | NewtonSqrt produces positive finite BL with indels |
| `test_optimize_method_newton_log_positive_with_indels` (3 cases)  | NewtonLog produces positive finite BL with indels  |

---

## Dispatch Zero Boundary

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_dispatch_zero_boundary.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_dispatch_zero_boundary.rs)

**Impl:**

- [`packages/treetime/src/commands/optimize/optimize_unified.rs`](../../packages/treetime/src/commands/optimize/optimize_unified.rs)
- [`packages/treetime/src/commands/optimize/method_newton.rs`](../../packages/treetime/src/commands/optimize/method_newton.rs)
- [`packages/treetime/src/commands/optimize/optimize_dense.rs`](../../packages/treetime/src/commands/optimize/optimize_dense.rs)
- [`packages/treetime/src/commands/optimize/partition_ops.rs`](../../packages/treetime/src/commands/optimize/partition_ops.rs)
- [`packages/treetime/src/commands/optimize/run.rs`](../../packages/treetime/src/commands/optimize/run.rs)

### End-to-end dispatch

| Test                                                                  | Purpose                                                      |
| --------------------------------------------------------------------- | ------------------------------------------------------------ |
| `test_dispatch_zero_boundary_k80_identical_sequences`                 | All 6 methods return t=0 on K80 identical sequences          |
| `test_dispatch_zero_boundary_non_unimodal_models_all_reach_zero`      | BrentSqrt reaches t=0 on every non-unimodal nucleotide model |
| `test_dispatch_zero_boundary_jc69_pre_dispatch_shortcut_reaches_zero` | `is_zero_branch_optimal` fires for JC69 identical sequences  |

### Reconcile helper contract

| Test                                                                           | Purpose                                                                       |
| ------------------------------------------------------------------------------ | ----------------------------------------------------------------------------- |
| `test_dispatch_zero_boundary_reconcile_positive_candidate_finds_positive_mode` | Positive candidate worse than zero -> grid returns the positive mode          |
| `test_dispatch_zero_boundary_reconcile_exact_zero_finds_positive_mode`         | Exact-zero candidate on multi-modal surface -> grid returns the positive mode |
| `test_dispatch_zero_boundary_reconcile_degenerate_site_passes_through`         | Degenerate site short-circuits via `all_sites_valid_at_zero` gate             |
| `test_dispatch_zero_boundary_reconcile_indel_count_positive_passes_through`    | `indel_count > 0` short-circuits via Poisson at zero                          |
| `test_dispatch_zero_boundary_reconcile_exact_zero_unimodal_passes_through`     | Unimodal model skips grid verification on exact-zero candidate                |
| `test_dispatch_zero_boundary_reconcile_exact_zero_indels_passes_through`       | Exact-zero candidate with indels skips grid verification                      |

### Downstream topology cleanup

| Test                                                                                 | Purpose                                                                                |
| ------------------------------------------------------------------------------------ | -------------------------------------------------------------------------------------- |
| `test_dispatch_zero_boundary_topology_cleanup_collects_k80_internal_edges` (4 cases) | `find_zero_optimal_internal_edges` collects both internal edges after K80 optimization |

### Inner-solver reproductions

| Test                                                                                 | Purpose                                                                             |
| ------------------------------------------------------------------------------------ | ----------------------------------------------------------------------------------- |
| `test_dispatch_zero_boundary_newton_inner_does_not_clamp_to_zero_on_dinh_matsen_k80` | Pin: `newton_inner` returns positive from 12 starting points on Dinh-Matsen surface |
| `test_dispatch_zero_boundary_newton_sqrt_inner_clamps_to_zero_on_dinh_matsen_k80`    | Reproduction: `newton_sqrt_inner` from t0=0.6 returns exactly 0 on the same surface |

---

## Optimization Method Step Clamping

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_optimize_method_step_clamping.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_optimize_method_step_clamping.rs)

**Impl:** [`packages/treetime/src/commands/optimize/method_newton.rs`](../../packages/treetime/src/commands/optimize/method_newton.rs)

| Test                                                                                | Purpose                                                 |
| ----------------------------------------------------------------------------------- | ------------------------------------------------------- |
| `test_optimize_method_step_clamping_sqrt_at_zero`                                   | At s=0 the bound equals -1.0 (coincides with t-space)   |
| `test_optimize_method_step_clamping_sqrt_analytical` (4 cases)                      | Known analytical values at representative s values      |
| `test_optimize_method_step_clamping_sqrt_always_negative` (7 cases)                 | Bound is always negative for non-negative s             |
| `test_optimize_method_step_clamping_sqrt_produces_unit_t_increase` (6 cases)        | Applying the bound step yields delta_t = 1.0            |
| `test_optimize_method_step_clamping_sqrt_large_s_asymptote`                         | For large s, bound approaches -1/(2s)                   |
| `test_optimize_method_step_clamping_sqrt_strictly_greater_than_minus_one` (4 cases) | For s > 0 the bound is strictly greater than -1.0       |
| `test_optimize_method_step_clamping_log_analytical` (4 cases)                       | Known analytical values at representative t values      |
| `test_optimize_method_step_clamping_log_always_negative` (6 cases)                  | Bound is always negative for positive t                 |
| `test_optimize_method_step_clamping_log_produces_unit_t_increase` (6 cases)         | Applying the bound step in u-space yields delta_t = 1.0 |
| `test_optimize_method_step_clamping_log_large_t_asymptote`                          | For large t, bound approaches -1/t                      |

---

## CLI Args

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_args.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_args.rs)

**Impl:** [`packages/treetime/src/commands/optimize/args.rs`](../../packages/treetime/src/commands/optimize/args.rs)

| Test                                               | Purpose                                                 |
| -------------------------------------------------- | ------------------------------------------------------- |
| `test_args_opt_method_kebab_case_parses` (6 cases) | Each BranchOptMethod kebab-case variant parses from CLI |
| `test_args_opt_method_default_is_brent_sqrt`       | Default --opt-method is BrentSqrt                       |
| `test_args_opt_method_rejects_unknown`             | Unknown --opt-method value rejected at parse time       |

---

## Zero Sequence Length

**Test:** [`packages/treetime/src/commands/optimize/__tests__/test_optimize_zero_sequence_length.rs`](../../packages/treetime/src/commands/optimize/__tests__/test_optimize_zero_sequence_length.rs)

**Impl:**

- [`packages/treetime/src/commands/optimize/optimize_unified.rs`](../../packages/treetime/src/commands/optimize/optimize_unified.rs)
- [`packages/treetime/src/commands/optimize/run.rs`](../../packages/treetime/src/commands/optimize/run.rs)

| Test                                                     | Purpose                                                      |
| -------------------------------------------------------- | ------------------------------------------------------------ |
| `test_optimize_zero_sequence_length_run_optimize_error`  | run_optimize_mixed returns error for zero-length partitions  |
| `test_optimize_zero_sequence_length_initial_guess_error` | initial_guess_mixed returns error for zero-length partitions |
