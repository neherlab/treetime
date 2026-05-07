# GTR Substitution Model Tests

[Back to index](README.md)

## Summary

| Category                                                              | Type            |
| --------------------------------------------------------------------- | --------------- |
| [Model construction (golden-master)](#golden-master-model-tests)      | Golden-master   |
| [Matrix exponentiation properties](#matrix-exponentiation-properties) | Property        |
| [Eigendecomposition properties](#eigendecomposition-properties)       | Property        |
| [Q matrix properties](#q-matrix-properties)                           | Property        |
| [Numerical stability](#numerical-stability-properties)                | Property        |
| [Branch length edge cases](#numerical-edge-cases)                     | Unit            |
| [Extreme parameter edge cases](#numerical-edge-cases)                 | Unit            |
| [Parameterized edge cases](#numerical-edge-cases)                     | Parameterized   |
| [Model hierarchy](#model-hierarchy-tests)                             | Unit            |
| [GTR output (write_gtr_json)](#gtr-output-tests)                      | Unit            |
| [Generator validation](#generator-validation)                         | Property        |
| [GTR inference (dense golden)](#gtr-inference-tests---dense)          | Golden-master   |
| [GTR inference (dense unit)](#gtr-inference-tests---dense)            | Unit            |
| [GTR inference (sparse)](#gtr-inference-tests---sparse)               | Unit            |
| [Inference contracts](#inference-contract-tests)                      | Parameterized   |
| [Dense-sparse cross-validation](#dense-sparse-cross-validation)       | Integration     |
| [Common functions](#common-functions)                                 | Unit            |
| [Site-specific GTR properties](#site-specific-gtr-tests)              | Property        |
| [Site-specific GTR unit+validation](#site-specific-gtr-tests)         | Unit            |
| [Site-specific GTR golden-master](#site-specific-gtr-tests)           | Golden-master   |
| [Site-specific GTR inference](#site-specific-gtr-inference-tests)     | Unit            |
| [Site-rate variation (inline)](#site-rate-variation)                  | Unit + Property |
| [Site-rate variation properties](#site-rate-variation-property-tests) | Property + Unit |
| [Jukes-Cantor distance correction](#jukes-cantor-distance-correction) | Mixed           |

Property tests use proptest with random inputs. Generator tests use smaller sample sizes.

---

## Golden-Master Model Tests

**Test:** [`packages/treetime/src/gtr/__tests__/test_gm_gtr.rs`](../../packages/treetime/src/gtr/__tests__/test_gm_gtr.rs)

**Impl:**

- [`packages/treetime/src/gtr/get_gtr.rs`](../../packages/treetime/src/gtr/get_gtr.rs)
- [`packages/treetime/src/gtr/gtr.rs`](../../packages/treetime/src/gtr/gtr.rs)

Validates Rust v1 GTR model construction against Python v0 reference outputs.

| Test                 | Model | Parameters                                                            |
| -------------------- | ----- | --------------------------------------------------------------------- |
| `test_gm_gtr_jc69`   | JC69  | default, mu_0_5, mu_2_0                                               |
| `test_gm_gtr_k80`    | K80   | default, kappa_0_5, kappa_1_0, kappa_2_0, kappa_5_0, mu_0_5_kappa_2_0 |
| `test_gm_gtr_f81`    | F81   | default, custom_pi, asymmetric_pi, mu_2_0                             |
| `test_gm_gtr_hky85`  | HKY85 | default, kappa_2_0, custom_pi, custom_all                             |
| `test_gm_gtr_t92`    | T92   | default, high_gc, low_gc, kappa_2_0, custom_all                       |
| `test_gm_gtr_tn93`   | TN93  | default, kappa1_0_5, kappa2_2_0, both_kappa, custom_pi, full_custom   |
| `test_gm_gtr_custom` | GTR   | uniform_rates, ts_tv_bias, asymmetric_pi, mu_scaling, random          |

**Tolerance:** 1e-14 (machine precision)
**Fixtures:** [`gm_gtr_inputs.json`](../../packages/treetime/src/gtr/__tests__/__fixtures__/gm_gtr_inputs.json), [`gm_gtr_outputs.json`](../../packages/treetime/src/gtr/__tests__/__fixtures__/gm_gtr_outputs.json)

---

## Matrix Exponentiation Properties

**Test:** [`packages/treetime/src/gtr/__tests__/test_prop_gtr_expqt.rs`](../../packages/treetime/src/gtr/__tests__/test_prop_gtr_expqt.rs)

**Impl:** [`packages/treetime/src/gtr/gtr.rs`](../../packages/treetime/src/gtr/gtr.rs)

Verifies P(t) = exp(Qt) mathematical invariants.

| Property                                            | Invariant                             |
| --------------------------------------------------- | ------------------------------------- |
| `test_prop_gtr_expqt_stochastic_columns`            | sum_i P\[i,j\]\(t\) = 1 for all j     |
| `test_prop_gtr_expqt_nonnegative`                   | P\[i,j\]\(t\) >= 0                    |
| `test_prop_gtr_expqt_bounded`                       | P\[i,j\]\(t\) <= 1                    |
| `test_prop_gtr_expqt_zero_is_identity`              | P\(0\) = I                            |
| `test_prop_gtr_expqt_equilibrium_limit`             | lim\(t->inf\) P\[i,j\]\(t\) = pi\[i\] |
| `test_prop_gtr_expqt_semigroup`                     | P(s+t) = P(s) \* P(t)                 |
| `test_prop_gtr_expqt_stationary_preserved`          | P(t) \* pi = pi                       |
| `test_prop_gtr_expqt_evolve_transpose_of_propagate` | Forward/backward transpose relation   |

---

## Eigendecomposition Properties

**Test:** [`packages/treetime/src/gtr/__tests__/test_prop_gtr_eigen.rs`](../../packages/treetime/src/gtr/__tests__/test_prop_gtr_eigen.rs)

**Impl:** [`packages/treetime/src/gtr/gtr.rs`](../../packages/treetime/src/gtr/gtr.rs)

Verifies eigendecomposition invariants.

| Property                                                           | Invariant              |
| ------------------------------------------------------------------ | ---------------------- |
| `test_prop_gtr_eigen_eigvals_nonpositive`                          | Re(lambda_i) <= 0      |
| `test_prop_gtr_eigen_eigvals_one_zero_eigenvalue`                  | Exactly one lambda = 0 |
| `test_prop_gtr_eigen_eigendecomposition_v_times_v_inv_is_identity` | V \* V^-1 = I          |

---

## Q Matrix Properties

**Test:** [`packages/treetime/src/gtr/__tests__/test_prop_gtr_q.rs`](../../packages/treetime/src/gtr/__tests__/test_prop_gtr_q.rs)

**Impl:**

- [`packages/treetime/src/gtr/gtr.rs`](../../packages/treetime/src/gtr/gtr.rs)
- [`packages/treetime/src/gtr/get_gtr.rs`](../../packages/treetime/src/gtr/get_gtr.rs)

Verifies rate matrix Q invariants.

| Property                              | Invariant                   |
| ------------------------------------- | --------------------------- |
| `test_prop_gtr_q_columns_sum_to_zero` | sum_i Q[i,j] = 0            |
| `test_prop_gtr_q_offdiag_nonnegative` | Q[i,j] >= 0 for i != j      |
| `test_prop_gtr_q_diag_nonpositive`    | Q[j,j] <= 0                 |
| `test_prop_gtr_q_detailed_balance`    | pi[j]*Q[i,j] = pi[i]*Q[j,i] |
| `test_prop_gtr_q_w_symmetric`         | W[i,j] = W[j,i]             |
| `test_prop_gtr_q_pi_sums_to_one`      | sum_i pi[i] = 1             |
| `test_prop_gtr_q_pi_positive`         | pi[i] > 0 for all i         |
| `test_prop_gtr_q_mu_scaling`          | exp(mu*Q*t) = exp(Q*(mu*t)) |

---

## Numerical Stability Properties

**Test:** [`packages/treetime/src/gtr/__tests__/test_prop_gtr_numerical.rs`](../../packages/treetime/src/gtr/__tests__/test_prop_gtr_numerical.rs)

**Impl:** [`packages/treetime/src/gtr/gtr.rs`](../../packages/treetime/src/gtr/gtr.rs)

Verifies no NaN/Inf and valid outputs.

| Property                                                 | Purpose                     |
| -------------------------------------------------------- | --------------------------- |
| `test_prop_gtr_numerical_no_nan_in_expqt`                | No NaN in exp(Qt)           |
| `test_prop_gtr_numerical_no_inf_in_expqt`                | No Inf in exp(Qt)           |
| `test_prop_gtr_numerical_propagate_profile_valid_output` | Non-negative, no NaN/Inf    |
| `test_prop_gtr_numerical_evolve_preserves_probability`   | Rows sum to 1, non-negative |

---

## Numerical Edge Cases

### Branch Length

**Test:** [`packages/treetime/src/gtr/__tests__/test_gtr_numerical_edge/test_gtr_numerical_edge_branch_length.rs`](../../packages/treetime/src/gtr/__tests__/test_gtr_numerical_edge/test_gtr_numerical_edge_branch_length.rs)

**Impl:**

- [`packages/treetime/src/gtr/gtr.rs`](../../packages/treetime/src/gtr/gtr.rs)
- [`packages/treetime/src/gtr/get_gtr.rs`](../../packages/treetime/src/gtr/get_gtr.rs)

| Test                                 | Scenario  | Expected                        |
| ------------------------------------ | --------- | ------------------------------- |
| `test_gtr_expqt_zero_branch`         | t = 0     | P(0) = I exactly                |
| `test_gtr_expqt_tiny_branch`         | t = 1e-10 | P ~ I, no underflow             |
| `test_gtr_expqt_small_branch_taylor` | t = 1e-6  | P ~ I + Qt (first-order)        |
| `test_gtr_expqt_large_branch`        | t = 100   | P ~ equilibrium (rows = pi)     |
| `test_gtr_expqt_very_large_branch`   | t = 1000  | Full equilibration, no overflow |

### Extreme Parameters

**Test:** [`packages/treetime/src/gtr/__tests__/test_gtr_numerical_edge/test_gtr_numerical_edge_extreme_parameters.rs`](../../packages/treetime/src/gtr/__tests__/test_gtr_numerical_edge/test_gtr_numerical_edge_extreme_parameters.rs)

**Impl:**

- [`packages/treetime/src/gtr/gtr.rs`](../../packages/treetime/src/gtr/gtr.rs)
- [`packages/treetime/src/gtr/get_gtr.rs`](../../packages/treetime/src/gtr/get_gtr.rs)

| Test                                    | Parameter | Edge Value                    | Notes            |
| --------------------------------------- | --------- | ----------------------------- | ---------------- |
| `test_gtr_k80_kappa_near_zero`          | kappa     | 0.01 (transitions suppressed) |                  |
| `test_gtr_k80_kappa_large`              | kappa     | 100.0 (transitions dominant)  |                  |
| `test_gtr_hky85_skewed_pi`              | pi        | [0.97, 0.01, 0.01, 0.01]      |                  |
| `test_gtr_hky85_nearly_uniform_pi`      | pi        | [0.24, 0.26, 0.25, 0.25]      |                  |
| `test_gtr_hky85_uniform_pi_matches_k80` | pi        | [0.25, 0.25, 0.25, 0.25]      | **1e-3** epsilon |
| `test_gtr_mu_very_small`                | mu        | 0.001 (slow evolution)        |                  |
| `test_gtr_mu_large`                     | mu        | 10.0 (fast evolution)         |                  |

### Parameterized Edge Cases

**Test:** [`packages/treetime/src/gtr/__tests__/test_gtr_numerical_edge/test_gtr_numerical_edge_parameterized.rs`](../../packages/treetime/src/gtr/__tests__/test_gtr_numerical_edge/test_gtr_numerical_edge_parameterized.rs)

**Impl:**

- [`packages/treetime/src/gtr/gtr.rs`](../../packages/treetime/src/gtr/gtr.rs)
- [`packages/treetime/src/gtr/get_gtr.rs`](../../packages/treetime/src/gtr/get_gtr.rs)

| Test                                 | Parameter Range         |
| ------------------------------------ | ----------------------- |
| `test_gtr_expqt_branch_length_range` | 0.0 to 1000.0           |
| `test_gtr_k80_kappa_range`           | 0.001 to 500.0          |
| `test_gtr_extreme_w_values`          | W entries 0.01 to 100.0 |

---

## Model Hierarchy Tests

**Test:** [`packages/treetime/src/gtr/__tests__/test_gtr_hierarchy/test_gtr_hierarchy_model_relationships.rs`](../../packages/treetime/src/gtr/__tests__/test_gtr_hierarchy/test_gtr_hierarchy_model_relationships.rs)

**Impl:**

- [`packages/treetime/src/gtr/get_gtr.rs`](../../packages/treetime/src/gtr/get_gtr.rs)
- [`packages/treetime/src/gtr/gtr.rs`](../../packages/treetime/src/gtr/gtr.rs)

Verifies nested model relationships:

```
JC69 ----> K80 ----> HKY85 ----> TN93 ----> GTR
  |                    ^
  +---> F81 ----------+
```

| Test                                           | Reduction    | Constraint               |
| ---------------------------------------------- | ------------ | ------------------------ |
| `test_gtr_jc69_equals_k80_kappa_1`             | JC69 = K80   | kappa = 1                |
| `test_gtr_jc69_equals_f81_uniform_pi`          | JC69 = F81   | pi = uniform             |
| `test_gtr_k80_equals_hky85_uniform_pi`         | K80 = HKY85  | pi = uniform             |
| `test_gtr_f81_equals_hky85_kappa_1`            | F81 = HKY85  | kappa = 1                |
| `test_gtr_hky85_equals_tn93_equal_transitions` | HKY85 = TN93 | kappa1 = 0.5, kappa2 = 1 |
| `test_gtr_tn93_equals_gtr_with_structured_w`   | TN93 = GTR   | W with TN93 structure    |

**Test:** [`packages/treetime/src/gtr/__tests__/test_gtr_hierarchy/test_gtr_hierarchy_additional_verification.rs`](../../packages/treetime/src/gtr/__tests__/test_gtr_hierarchy/test_gtr_hierarchy_additional_verification.rs)

**Impl:**

- [`packages/treetime/src/gtr/get_gtr.rs`](../../packages/treetime/src/gtr/get_gtr.rs)
- [`packages/treetime/src/gtr/gtr.rs`](../../packages/treetime/src/gtr/gtr.rs)

| Test                                  | Verification                     |
| ------------------------------------- | -------------------------------- |
| `test_gtr_jc69_q_symmetric`           | JC69 has symmetric Q             |
| `test_gtr_k80_kappa_1_q_symmetric`    | K80(kappa=1) has symmetric Q     |
| `test_gtr_mu_does_not_affect_q_shape` | Q normalized independently of mu |

---

## Generator Validation

**Test:** [`packages/treetime/src/gtr/__tests__/generators.rs`](../../packages/treetime/src/gtr/__tests__/generators.rs)

**Impl:** [`packages/treetime/src/gtr/gtr.rs`](../../packages/treetime/src/gtr/gtr.rs)

Validates that proptest generators produce valid outputs.

| Property                                             | Invariant                         |
| ---------------------------------------------------- | --------------------------------- |
| `test_prop_generators_arb_pi_nuc_sums_to_one`        | pi sums to 1, all positive        |
| `test_prop_generators_arb_pi_aa_sums_to_one`         | pi sums to 1, all positive        |
| `test_prop_generators_arb_w_nuc_symmetric_and_valid` | W symmetric, zero diagonal        |
| `test_prop_generators_arb_w_aa_symmetric_and_valid`  | W symmetric, zero diagonal        |
| `test_prop_generators_arb_branch_len_range`          | Branch length in valid range      |
| `test_prop_generators_arb_profile_nuc_valid`         | Profile rows sum to 1             |
| `test_prop_generators_arb_gtr_nuc_valid`             | GTR model passes basic invariants |

---

## GTR Inference Tests - Dense

### Golden-Master

**Test:** [`packages/treetime/src/gtr/infer_gtr/__tests__/test_gm_infer_gtr_dense.rs`](../../packages/treetime/src/gtr/infer_gtr/__tests__/test_gm_infer_gtr_dense.rs)

**Impl:** [`packages/treetime/src/gtr/infer_gtr/dense.rs`](../../packages/treetime/src/gtr/infer_gtr/dense.rs)

| Test                                | Datasets                                                                                                   |
| ----------------------------------- | ---------------------------------------------------------------------------------------------------------- |
| `test_gm_infer_gtr_dense_synthetic` | simple_4taxa, star_topology, caterpillar, multiple_mutations, varying_branch_lengths, large_branchy_uneven |
| `test_gm_infer_gtr_dense_real`      | dengue_20, ebola_20, flu_h3n2_20, rsv_a_20, lassa_L_50, mpox_clade_ii_20, tb_20                            |

**Tolerances:**

- Synthetic: 1e-8 (short sequences)
- Real: 1e-6 (BLAS drift scales with sequence length)

### Unit Tests

**Test:** [`packages/treetime/src/gtr/infer_gtr/__tests__/test_dense.rs`](../../packages/treetime/src/gtr/infer_gtr/__tests__/test_dense.rs)

**Impl:**

- [`packages/treetime/src/gtr/infer_gtr/dense.rs`](../../packages/treetime/src/gtr/infer_gtr/dense.rs)
- [`packages/treetime/src/gtr/infer_gtr/common.rs`](../../packages/treetime/src/gtr/infer_gtr/common.rs)

| Test                                 | Purpose                                     |
| ------------------------------------ | ------------------------------------------- |
| `test_uniform_sequences`             | Identical sequences: off-diagonal nij small |
| `test_single_mutation`               | One A->C mutation: nij carries signal       |
| `test_zero_branch_lengths`           | BL = 0 (clamped): Ti small, bounded nij     |
| `test_zero_branch_lengths_unclamped` | BL = 0, expQt = I: Ti = 0, nij diagonal     |
| `test_produces_valid_model`          | W symmetric, pi normalized, mu > 0          |

---

## GTR Inference Tests - Sparse

**Test:** [`packages/treetime/src/gtr/infer_gtr/__tests__/test_sparse.rs`](../../packages/treetime/src/gtr/infer_gtr/__tests__/test_sparse.rs)

**Impl:**

- [`packages/treetime/src/gtr/infer_gtr/sparse.rs`](../../packages/treetime/src/gtr/infer_gtr/sparse.rs)
- [`packages/treetime/src/gtr/infer_gtr/common.rs`](../../packages/treetime/src/gtr/infer_gtr/common.rs)

| Test                              | Purpose                                      |
| --------------------------------- | -------------------------------------------- |
| `test_get_mutation_counts_sparse` | Integer mutation counts from Fitch parsimony |
| `test_infer_gtr_sparse`           | GTR inference from sparse representation     |

---

## Inference Contract Tests

**Test:** [`packages/treetime/src/gtr/infer_gtr/__tests__/test_contract.rs`](../../packages/treetime/src/gtr/infer_gtr/__tests__/test_contract.rs)

**Impl:**

- [`packages/treetime/src/gtr/infer_gtr/dense.rs`](../../packages/treetime/src/gtr/infer_gtr/dense.rs)
- [`packages/treetime/src/gtr/infer_gtr/sparse.rs`](../../packages/treetime/src/gtr/infer_gtr/sparse.rs)

Verifies mutation count invariants across dense and sparse paths.

| Test                                                  | Purpose                                  | Notes                                   |
| ----------------------------------------------------- | ---------------------------------------- | --------------------------------------- |
| `test_nij_orientation_dense`                          | A->C: nij[C,A] >> nij[A,C] (dense)       |                                         |
| `test_nij_orientation_sparse`                         | A->C: nij[C,A] = 1, nij[A,C] = 0         |                                         |
| `test_ti_scaling_sparse` (2 cases)                    | Ti proportional to branch length         |                                         |
| `test_ti_scaling_dense` (2 cases)                     | Ti proportional to branch length         |                                         |
| `test_dense_sparse_consistency`                       | Dense ~ sparse for unambiguous sequences | **1e-1** nij diff, **1e-2** ti rel diff |
| `test_root_state_dense`                               | Root reflects ancestral composition      |                                         |
| `test_root_state_sparse`                              | Root reflects Fitch consensus            |                                         |
| `test_nij_orientation_multiple_mutations_sparse`      | Multiple mutation directions correct     |                                         |
| `test_nij_accumulation_dense`                         | nij accumulates across edges             |                                         |
| `test_ti_proportional_to_composition_sparse`          | Ti proportional to state frequency       |                                         |
| `test_dense_sparse_nij_direction_agreement`           | Both agree on dominant mutation cell     |                                         |
| `test_root_state_total_equals_alignment_length_dense` | sum(root_state) = seq_length             |                                         |
| `test_nij_diagonal_zero`                              | nij[i,i] = 0 (dense and sparse)          |                                         |
| `test_nij_non_negative`                               | nij[i,j] >= 0 (dense and sparse)         |                                         |
| `test_ti_non_negative`                                | Ti[i] >= 0 (dense and sparse)            |                                         |

---

## Dense-Sparse Cross-Validation

**Test:** [`packages/treetime/src/gtr/infer_gtr/__tests__/test_contract_dense_sparse_real.rs`](../../packages/treetime/src/gtr/infer_gtr/__tests__/test_contract_dense_sparse_real.rs)

**Impl:**

- [`packages/treetime/src/gtr/infer_gtr/dense.rs`](../../packages/treetime/src/gtr/infer_gtr/dense.rs)
- [`packages/treetime/src/gtr/infer_gtr/sparse.rs`](../../packages/treetime/src/gtr/infer_gtr/sparse.rs)

| Test                                  | Datasets                                                                        |
| ------------------------------------- | ------------------------------------------------------------------------------- |
| `test_contract_dense_sparse_real_gtr` | flu_h3n2_20, ebola_20, rsv_a_20, dengue_20, tb_20, lassa_L_50, mpox_clade_ii_20 |

**Thresholds:**

- pi cosine > 0.997
- W relative Frobenius < 0.21
- mu relative diff < 0.23

---

## Common Functions

**Test:** [`packages/treetime/src/gtr/infer_gtr/__tests__/test_common.rs`](../../packages/treetime/src/gtr/infer_gtr/__tests__/test_common.rs)

**Impl:** [`packages/treetime/src/gtr/infer_gtr/common.rs`](../../packages/treetime/src/gtr/infer_gtr/common.rs)

| Test                              | Purpose                                  |
| --------------------------------- | ---------------------------------------- |
| `test_infer_gtr_impl`             | W, pi, mu inference from mutation counts |
| `test_distance_zero_to_uniform`   | Distance: zero -> uniform = 0.5          |
| `test_distance_uniform_to_skewed` | Distance: uniform -> skewed              |
| `test_distance_small_difference`  | Distance: small perturbations            |
| `test_distance_tiny_difference`   | Distance: tiny perturbations             |
| `test_distance_identical`         | Distance: identity = 0                   |

---

## GTR Output Tests

**Test:** [`packages/treetime/src/gtr/__tests__/test_write_gtr_json.rs`](../../packages/treetime/src/gtr/__tests__/test_write_gtr_json.rs)

**Impl:** [`packages/treetime/src/gtr/get_gtr.rs`](../../packages/treetime/src/gtr/get_gtr.rs)

| Test                                               | Purpose                                                 |
| -------------------------------------------------- | ------------------------------------------------------- |
| `test_write_gtr_json_filename` (3 cases)           | Qualifier produces correct filename                     |
| `test_write_gtr_json_both_partitions_no_overwrite` | Both sparse and dense qualifiers produce separate files |

---

## Jukes-Cantor Distance Correction

**Test:** [`packages/treetime/src/gtr/jc_distance.rs`](../../packages/treetime/src/gtr/jc_distance.rs)

**Impl:** [`packages/treetime/src/gtr/jc_distance.rs`](../../packages/treetime/src/gtr/jc_distance.rs)

| Test                                                           | Purpose                                                      |
| -------------------------------------------------------------- | ------------------------------------------------------------ |
| `test_jukes_cantor_distance_zero_p_returns_zero`               | `p = 0` returns exactly `0` for any alphabet                 |
| `test_jukes_cantor_distance_negative_p_clamped_to_zero`        | Negative `p` defensively clamps to 0                         |
| `test_jukes_cantor_distance_known_values` (6 cases)            | Analytical values for k in {4, 20} at representative p       |
| `test_jukes_cantor_distance_saturation_is_finite` (5 cases)    | Saturation and beyond produce finite positive distances      |
| `test_jukes_cantor_distance_saturation_cap_order_of_magnitude` | Saturation cap lies in expected bracket                      |
| `test_jukes_cantor_distance_always_at_least_p`                 | Property: d(p) >= p across the valid range                   |
| `test_jukes_cantor_distance_monotonic_in_p`                    | Property: d is non-decreasing in p                           |
| `test_jukes_cantor_distance_small_p_approaches_p`              | Taylor behaviour: d -> p as p -> 0                           |
| `test_jukes_cantor_distance_issue_documented_error`            | Correction sizes match the 7% / 22% figures at p = 0.1, 0.25 |

---

## Site-Specific GTR Tests

### Property Tests

**Test:** [`packages/treetime/src/gtr/__tests__/test_prop_gtr_site_specific.rs`](../../packages/treetime/src/gtr/__tests__/test_prop_gtr_site_specific.rs)

**Impl:** [`packages/treetime/src/gtr/gtr_site_specific.rs`](../../packages/treetime/src/gtr/gtr_site_specific.rs)

| Test                                                        | Purpose                                              | Notes              |
| ----------------------------------------------------------- | ---------------------------------------------------- | ------------------ |
| `test_prop_gtr_site_specific_expqt_column_stochastic`       | P_a(t) columns sum to 1 per site                     |                    |
| `test_prop_gtr_site_specific_expqt_identity_at_zero`        | P_a(0) = I per site                                  |                    |
| `test_prop_gtr_site_specific_expqt_nonnegative`             | P_a(t) entries non-negative                          |                    |
| `test_prop_gtr_site_specific_expqt_semigroup`               | P_a(s+t) = P_a(s) \* P_a(t) per site                 |                    |
| `test_prop_gtr_site_specific_expqt_convergence`             | P_a(t) -> pi_a as t -> infinity                      |                    |
| `test_prop_gtr_site_specific_propagate_profile_valid`       | propagate_profile output finite and non-negative     |                    |
| `test_prop_gtr_site_specific_evolve_valid`                  | evolve output finite and non-negative                |                    |
| `test_prop_gtr_site_specific_equilibrium_fixed_point`       | evolve(pi, t) = pi (equilibrium is fixed point)      |                    |
| `test_prop_gtr_site_specific_interpolation_accuracy`        | Interpolated expQt within 1e-2 of direct computation | **1e-2** tolerance |
| `test_prop_gtr_site_specific_average_rate_positive`         | Per-site average rate is positive                    |                    |
| `test_prop_gtr_site_specific_expqt_bounded`                 | P_a(t) entries bounded by 1                          |                    |
| `test_prop_gtr_site_specific_stationary_preserved`          | P_a(t) @ pi_a = pi_a (right eigenvector)             |                    |
| `test_prop_gtr_site_specific_evolve_transpose_of_propagate` | evolve = profile @ P^T, propagate = profile @ P      |                    |
| `test_prop_gtr_site_specific_no_nan`                        | No NaN in expQt                                      |                    |
| `test_prop_gtr_site_specific_no_inf`                        | No Inf in expQt                                      |                    |
| `test_prop_gtr_site_specific_evolve_preserves_probability`  | evolve preserves row sums                            |                    |
| `test_prop_gtr_site_specific_approx_column_stochastic`      | Interpolation preserves column stochasticity         |                    |
| `test_prop_gtr_site_specific_approx_nonnegative`            | Interpolation preserves non-negativity               |                    |
| `test_prop_gtr_site_specific_approx_equilibrium`            | Interpolation preserves equilibrium                  |                    |

### Unit + Validation Tests

**Test:** [`packages/treetime/src/gtr/__tests__/test_prop_gtr_site_specific.rs`](../../packages/treetime/src/gtr/__tests__/test_prop_gtr_site_specific.rs) (same file)

**Impl:** [`packages/treetime/src/gtr/gtr_site_specific.rs`](../../packages/treetime/src/gtr/gtr_site_specific.rs)

| Test                                                | Purpose                                       |
| --------------------------------------------------- | --------------------------------------------- |
| `test_gtr_site_specific_different_sites_differ`     | Different pi produces different P(t) per site |
| `test_gtr_site_specific_uniform_matches_standard`   | Uniform-pi site-specific matches standard GTR |
| `test_gtr_site_specific_rejects_zero_pi_column`     | Constructor rejects zero-sum pi column        |
| `test_gtr_site_specific_rejects_negative_mu`        | Constructor rejects negative mu               |
| `test_gtr_site_specific_rejects_dimension_mismatch` | Constructor rejects mu length mismatch        |

### Golden-Master Tests

**Test:** [`packages/treetime/src/gtr/__tests__/test_gm_gtr_site_specific.rs`](../../packages/treetime/src/gtr/__tests__/test_gm_gtr_site_specific.rs)

**Impl:**

- [`packages/treetime/src/gtr/gtr_site_specific.rs`](../../packages/treetime/src/gtr/gtr_site_specific.rs)
- [`packages/treetime/src/gtr/infer_gtr/site_specific.rs`](../../packages/treetime/src/gtr/infer_gtr/site_specific.rs)

| Test                                         | Tolerance | Purpose                                            | Notes                         |
| -------------------------------------------- | --------- | -------------------------------------------------- | ----------------------------- |
| `test_gm_gtr_site_specific_expqt`            | 1e-10     | Eigenvalues and expQt against v0 oracle            |                               |
| `test_gm_gtr_site_specific_propagate_evolve` | 1e-10     | propagate_profile and evolve against v0 oracle     |                               |
| `test_gm_gtr_site_specific_infer`            | **1e-3**  | Inference pi and W ratios against v0 oracle        | Circular inference distortion |
| `test_gm_gtr_site_specific_approximate`      | 1e-8      | Interpolated expQt against v0 interpolation oracle |                               |

---

## Site-Specific GTR Inference Tests

**Test:** [`packages/treetime/src/gtr/infer_gtr/__tests__/test_site_specific.rs`](../../packages/treetime/src/gtr/infer_gtr/__tests__/test_site_specific.rs)

**Impl:** [`packages/treetime/src/gtr/infer_gtr/site_specific.rs`](../../packages/treetime/src/gtr/infer_gtr/site_specific.rs)

| Test                                                | Purpose                                         | Notes                          |
| --------------------------------------------------- | ----------------------------------------------- | ------------------------------ |
| `test_infer_gtr_site_specific_recovers_parameters`  | Inference from synthetic data recovers W and pi | **1e-3** pi, **1e-2** W ratios |
| `test_infer_gtr_site_specific_produces_valid_model` | Inferred model produces valid P(t) matrices     |                                |

---

## Site-Rate Variation

**Test:** [`packages/treetime/src/gtr/site_rate_variation.rs`](../../packages/treetime/src/gtr/site_rate_variation.rs) (inline `#[cfg(test)]`)

**Impl:** [`packages/treetime/src/gtr/site_rate_variation.rs`](../../packages/treetime/src/gtr/site_rate_variation.rs)

| Test                                                      | Purpose                                                        | Notes               |
| --------------------------------------------------------- | -------------------------------------------------------------- | ------------------- |
| `test_discrete_gamma_rates_mean_one` (7 cases)            | Discrete gamma rates have mean 1.0 across alpha/K combinations |                     |
| `test_discrete_gamma_rates_single_category`               | Single category returns [1.0]                                  |                     |
| `test_discrete_gamma_rates_sorted_ascending`              | Rates are sorted ascending                                     |                     |
| `test_discrete_gamma_rates_all_positive`                  | All rates positive                                             |                     |
| `test_discrete_gamma_rates_high_alpha_approaches_uniform` | Large alpha produces near-uniform rates                        |                     |
| `test_discrete_gamma_rates_low_alpha_wide_spread`         | Small alpha produces wide rate spread (ratio > 50)             |                     |
| `test_discrete_gamma_rates_invalid_alpha`                 | Invalid alpha values rejected with error                       |                     |
| `test_discrete_gamma_rates_invalid_categories`            | Zero categories rejected                                       |                     |
| `test_discrete_gamma_rates_reference_alpha_1_k4`          | Reference values for alpha=1 K=4 match Yang 1994               |                     |
| `test_prop_discrete_gamma_rates_mean_one`                 | Property: mean is 1.0 across random alpha/K                    | **ignored** (flaky) |
| `test_prop_discrete_gamma_rates_positive_sorted`          | Property: rates positive and sorted                            | **ignored** (flaky) |

---

## Site-Rate Variation Property Tests

**Test:** [`packages/treetime/src/gtr/__tests__/test_prop_gtr_site_rates.rs`](../../packages/treetime/src/gtr/__tests__/test_prop_gtr_site_rates.rs)

**Impl:** [`packages/treetime/src/gtr/gtr.rs`](../../packages/treetime/src/gtr/gtr.rs)

| Test                                              | Purpose                                                             |
| ------------------------------------------------- | ------------------------------------------------------------------- |
| `test_prop_expqt_with_rate_one_equals_expqt`      | expQt_with_rate(t, 1.0) equals expQt(t)                             |
| `test_prop_expqt_with_rate_column_stochastic`     | expQt_with_rate produces column-stochastic matrices                 |
| `test_prop_expqt_with_rate_non_negative`          | expQt_with_rate entries are non-negative                            |
| `test_prop_expqt_with_rate_semigroup`             | P(r*s) * P(r*t) = P(r*(s+t)) semigroup property                     |
| `test_prop_propagate_uniform_rates_equals_scalar` | propagate_profile with uniform site_rates equals scalar propagation |
| `test_prop_evolve_uniform_rates_equals_scalar`    | evolve with uniform site_rates equals scalar evolution              |
| `test_prop_propagate_per_site_matches_individual` | Per-site propagate matches row-wise individual expQt_with_rate      |
| `test_prop_evolve_per_site_matches_individual`    | Per-site evolve matches row-wise individual expQt_with_rate         |
| `test_prop_propagate_per_site_identity_at_zero`   | Per-site propagate at t=0 is identity                               |
| `test_prop_evolve_per_site_equilibrium`           | Per-site evolve at large t converges to equilibrium (pi)            |
| `test_higher_rate_more_divergence`                | Higher site rate produces more divergence                           |
| `test_rate_zero_is_identity`                      | Rate 0 produces identity matrix                                     |
| `test_site_rates_lifecycle`                       | set_site_rates / has_site_rates / clear_site_rates lifecycle        |

---

## Shared Test Support

**Test:** [`packages/treetime/src/gtr/__tests__/prop_support.rs`](../../packages/treetime/src/gtr/__tests__/prop_support.rs) - Proptest assertion helpers (`prop_assert_columns_sum_to`, `prop_assert_rows_sum_to`, `prop_assert_detailed_balance`)

**Test:** [`packages/treetime/src/gtr/__tests__/test_gtr_numerical_edge/test_gtr_numerical_edge_support.rs`](../../packages/treetime/src/gtr/__tests__/test_gtr_numerical_edge/test_gtr_numerical_edge_support.rs) - Stochastic matrix assertion helper (`assert_stochastic_matrix`)

**Test:** [`packages/treetime/src/gtr/__tests__/site_specific_support.rs`](../../packages/treetime/src/gtr/__tests__/site_specific_support.rs) - `simulate_counts()`, `value_to_array2()`, `value_to_array3()`
