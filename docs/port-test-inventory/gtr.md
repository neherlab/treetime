# GTR Substitution Model Tests

[Back to index](_index.md)

## Summary

| Category                           | Files  | Tests     | Type          |
| ---------------------------------- | ------ | --------- | ------------- |
| Model construction (golden-master) | 1      | 33        | Golden-master |
| Matrix exponentiation properties   | 1      | 8 x 256   | Property      |
| Eigendecomposition properties      | 1      | 3 x 256   | Property      |
| Q matrix properties                | 1      | 8 x 256   | Property      |
| Numerical stability                | 1      | 4 x 256   | Property      |
| Branch length edge cases           | 1      | 5         | Unit          |
| Extreme parameter edge cases       | 1      | 7         | Unit          |
| Parameterized edge cases           | 1      | 3         | Unit + rstest |
| Model hierarchy                    | 2      | 9         | Unit          |
| Generator validation               | 1      | 7 x 64    | Property      |
| GTR inference (dense golden)       | 1      | 13        | Golden-master |
| GTR inference (dense unit)         | 1      | 5         | Unit          |
| GTR inference (sparse)             | 1      | 2         | Unit          |
| Inference contracts                | 1      | 15        | Unit + rstest |
| Dense-sparse cross-validation      | 1      | 7         | Integration   |
| Common functions                   | 1      | 6         | Unit          |
| **Total**                          | **16** | **~6400** | Mixed         |

---

## Golden-Master Model Tests

**File:** [`test_gm_gtr.rs`](../../packages/treetime/src/gtr/__tests__/test_gm_gtr.rs)

Validates Rust v1 GTR model construction against Python v0 reference outputs.

| Test                 | Model | Cases | Parameters                                                            |
| -------------------- | ----- | ----- | --------------------------------------------------------------------- |
| `test_gm_gtr_jc69`   | JC69  | 3     | default, mu_0_5, mu_2_0                                               |
| `test_gm_gtr_k80`    | K80   | 6     | default, kappa_0_5, kappa_1_0, kappa_2_0, kappa_5_0, mu_0_5_kappa_2_0 |
| `test_gm_gtr_f81`    | F81   | 4     | default, custom_pi, asymmetric_pi, mu_2_0                             |
| `test_gm_gtr_hky85`  | HKY85 | 4     | default, kappa_2_0, custom_pi, custom_all                             |
| `test_gm_gtr_t92`    | T92   | 5     | default, high_gc, low_gc, kappa_2_0, custom_all                       |
| `test_gm_gtr_tn93`   | TN93  | 6     | default, kappa1_0_5, kappa2_2_0, both_kappa, custom_pi, full_custom   |
| `test_gm_gtr_custom` | GTR   | 5     | uniform_rates, ts_tv_bias, asymmetric_pi, mu_scaling, random          |

**Tolerance:** 1e-14 (machine precision)
**Fixtures:** `__fixtures__/gm_gtr_inputs.json`, `__fixtures__/gm_gtr_outputs.json`

---

## Matrix Exponentiation Properties

**File:** [`test_prop_gtr_expqt.rs`](../../packages/treetime/src/gtr/__tests__/test_prop_gtr_expqt.rs)

Verifies P(t) = exp(Qt) mathematical invariants (256 cases each).

| Property                                            | Invariant                           |
| --------------------------------------------------- | ----------------------------------- |
| `test_prop_gtr_expqt_stochastic_columns`            | sum_i P[i,j](t) = 1 for all j       |
| `test_prop_gtr_expqt_nonnegative`                   | P[i,j](t) >= 0                      |
| `test_prop_gtr_expqt_bounded`                       | P[i,j](t) <= 1                      |
| `test_prop_gtr_expqt_zero_is_identity`              | P(0) = I                            |
| `test_prop_gtr_expqt_equilibrium_limit`             | lim(t->inf) P[i,j](t) = pi[i]       |
| `test_prop_gtr_expqt_semigroup`                     | P(s+t) = P(s) \* P(t)               |
| `test_prop_gtr_expqt_stationary_preserved`          | P(t) \* pi = pi                     |
| `test_prop_gtr_expqt_evolve_transpose_of_propagate` | Forward/backward transpose relation |

---

## Eigendecomposition Properties

**File:** [`test_prop_gtr_eigen.rs`](../../packages/treetime/src/gtr/__tests__/test_prop_gtr_eigen.rs)

Verifies eigendecomposition invariants (256 cases each).

| Property                                                           | Invariant              |
| ------------------------------------------------------------------ | ---------------------- |
| `test_prop_gtr_eigen_eigvals_nonpositive`                          | Re(lambda_i) <= 0      |
| `test_prop_gtr_eigen_eigvals_one_zero_eigenvalue`                  | Exactly one lambda = 0 |
| `test_prop_gtr_eigen_eigendecomposition_v_times_v_inv_is_identity` | V \* V^-1 = I          |

---

## Q Matrix Properties

**File:** [`test_prop_gtr_q.rs`](../../packages/treetime/src/gtr/__tests__/test_prop_gtr_q.rs)

Verifies rate matrix Q invariants (256 cases each).

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

**File:** [`test_prop_gtr_numerical.rs`](../../packages/treetime/src/gtr/__tests__/test_prop_gtr_numerical.rs)

Verifies no NaN/Inf and valid outputs (256 cases each).

| Property                                                 | Purpose                     |
| -------------------------------------------------------- | --------------------------- |
| `test_prop_gtr_numerical_no_nan_in_expqt`                | No NaN in exp(Qt)           |
| `test_prop_gtr_numerical_no_inf_in_expqt`                | No Inf in exp(Qt)           |
| `test_prop_gtr_numerical_propagate_profile_valid_output` | Non-negative, no NaN/Inf    |
| `test_prop_gtr_numerical_evolve_preserves_probability`   | Rows sum to 1, non-negative |

---

## Numerical Edge Cases

### Branch Length

**File:** [`test_gtr_numerical_edge_branch_length.rs`](../../packages/treetime/src/gtr/__tests__/test_gtr_numerical_edge/test_gtr_numerical_edge_branch_length.rs)

| Test                                 | Scenario  | Expected                        |
| ------------------------------------ | --------- | ------------------------------- |
| `test_gtr_expqt_zero_branch`         | t = 0     | P(0) = I exactly                |
| `test_gtr_expqt_tiny_branch`         | t = 1e-10 | P ~ I, no underflow             |
| `test_gtr_expqt_small_branch_taylor` | t = 1e-6  | P ~ I + Qt (first-order)        |
| `test_gtr_expqt_large_branch`        | t = 100   | P ~ equilibrium (rows = pi)     |
| `test_gtr_expqt_very_large_branch`   | t = 1000  | Full equilibration, no overflow |

### Extreme Parameters

**File:** [`test_gtr_numerical_edge_extreme_parameters.rs`](../../packages/treetime/src/gtr/__tests__/test_gtr_numerical_edge/test_gtr_numerical_edge_extreme_parameters.rs)

| Test                                    | Parameter | Edge Value                    |
| --------------------------------------- | --------- | ----------------------------- |
| `test_gtr_k80_kappa_near_zero`          | kappa     | 0.01 (transitions suppressed) |
| `test_gtr_k80_kappa_large`              | kappa     | 100.0 (transitions dominant)  |
| `test_gtr_hky85_skewed_pi`              | pi        | [0.97, 0.01, 0.01, 0.01]      |
| `test_gtr_hky85_nearly_uniform_pi`      | pi        | [0.24, 0.26, 0.25, 0.25]      |
| `test_gtr_hky85_uniform_pi_matches_k80` | pi        | [0.25, 0.25, 0.25, 0.25]      |
| `test_gtr_mu_very_small`                | mu        | 0.001 (slow evolution)        |
| `test_gtr_mu_large`                     | mu        | 10.0 (fast evolution)         |

### Parameterized Edge Cases

**File:** [`test_gtr_numerical_edge_parameterized.rs`](../../packages/treetime/src/gtr/__tests__/test_gtr_numerical_edge/test_gtr_numerical_edge_parameterized.rs)

| Test                                 | Parameter Range         | Cases |
| ------------------------------------ | ----------------------- | ----- |
| `test_gtr_expqt_branch_length_range` | 0.0 to 1000.0           | 7     |
| `test_gtr_k80_kappa_range`           | 0.001 to 500.0          | 6     |
| `test_gtr_extreme_w_values`          | W entries 0.01 to 100.0 | 1     |

---

## Model Hierarchy Tests

**File:** [`test_gtr_hierarchy_model_relationships.rs`](../../packages/treetime/src/gtr/__tests__/test_gtr_hierarchy/test_gtr_hierarchy_model_relationships.rs)

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

**File:** [`test_gtr_hierarchy_additional_verification.rs`](../../packages/treetime/src/gtr/__tests__/test_gtr_hierarchy/test_gtr_hierarchy_additional_verification.rs)

| Test                                  | Verification                     |
| ------------------------------------- | -------------------------------- |
| `test_gtr_jc69_q_symmetric`           | JC69 has symmetric Q             |
| `test_gtr_k80_kappa_1_q_symmetric`    | K80(kappa=1) has symmetric Q     |
| `test_gtr_mu_does_not_affect_q_shape` | Q normalized independently of mu |

---

## Generator Validation

**File:** [`generators.rs`](../../packages/treetime/src/gtr/__tests__/generators.rs)

Validates that proptest generators produce valid outputs (64 cases each).

| Property                                             | Invariant                         |
| ---------------------------------------------------- | --------------------------------- |
| `test_prop_generators_arb_pi_nuc_sums_to_one`        | pi sums to 1, all positive        |
| `test_prop_generators_arb_pi_aa_sums_to_one`         | pi sums to 1, all positive        |
| `test_prop_generators_arb_w_nuc_symmetric_and_valid` | W symmetric, zero diagonal        |
| `test_prop_generators_arb_w_aa_symmetric_and_valid`  | W symmetric, zero diagonal        |
| `test_prop_generators_arb_branch_len_range`          | Branch length in valid range      |
| `test_prop_generators_arb_profile_nuc_valid`         | Profile rows sum to 1             |
| `test_prop_generators_arb_gtr_nuc_valid`             | GTR model passes basic invariants |

Generators (not tests):

| Generator         | Generates                         |
| ----------------- | --------------------------------- |
| `arb_pi_nuc`      | Nucleotide equilibrium freqs      |
| `arb_pi_aa`       | Amino acid equilibrium freqs      |
| `arb_w_nuc`       | Nucleotide exchangeability matrix |
| `arb_w_aa`        | Amino acid exchangeability matrix |
| `arb_branch_len`  | Log-uniform branch lengths        |
| `arb_profile_nuc` | Probability profiles              |
| `arb_gtr_nuc`     | Valid nucleotide GTR models       |

---

## GTR Inference Tests - Dense

### Golden-Master

**File:** [`test_gm_infer_gtr_dense.rs`](../../packages/treetime/src/gtr/infer_gtr/__tests__/test_gm_infer_gtr_dense.rs)

| Test                                | Cases | Datasets                                                                                                   |
| ----------------------------------- | ----- | ---------------------------------------------------------------------------------------------------------- |
| `test_gm_infer_gtr_dense_synthetic` | 6     | simple_4taxa, star_topology, caterpillar, multiple_mutations, varying_branch_lengths, large_branchy_uneven |
| `test_gm_infer_gtr_dense_real`      | 7     | dengue_20, ebola_20, flu_h3n2_20, rsv_a_20, lassa_L_50, mpox_clade_ii_20, tb_20                            |

**Tolerances:**

- Synthetic: 1e-8 (short sequences)
- Real: 1e-6 (BLAS drift scales with sequence length)

### Unit Tests

**File:** [`test_dense.rs`](../../packages/treetime/src/gtr/infer_gtr/__tests__/test_dense.rs)

| Test                                 | Purpose                                     |
| ------------------------------------ | ------------------------------------------- |
| `test_uniform_sequences`             | Identical sequences: off-diagonal nij small |
| `test_single_mutation`               | One A->C mutation: nij carries signal       |
| `test_zero_branch_lengths`           | BL = 0 (clamped): Ti small, bounded nij     |
| `test_zero_branch_lengths_unclamped` | BL = 0, expQt = I: Ti = 0, nij diagonal     |
| `test_produces_valid_model`          | W symmetric, pi normalized, mu > 0          |

---

## GTR Inference Tests - Sparse

**File:** [`test_sparse.rs`](../../packages/treetime/src/gtr/infer_gtr/__tests__/test_sparse.rs)

| Test                              | Purpose                                      |
| --------------------------------- | -------------------------------------------- |
| `test_get_mutation_counts_sparse` | Integer mutation counts from Fitch parsimony |
| `test_infer_gtr_sparse`           | GTR inference from sparse representation     |

---

## Inference Contract Tests

**File:** [`test_contract.rs`](../../packages/treetime/src/gtr/infer_gtr/__tests__/test_contract.rs)

Verifies mutation count invariants across dense and sparse paths.

| Test                                                  | Purpose                                    |
| ----------------------------------------------------- | ------------------------------------------ |
| `test_nij_orientation_dense`                          | A->C: nij[C,A] >> nij[A,C] (dense)         |
| `test_nij_orientation_sparse`                         | A->C: nij[C,A] = 1, nij[A,C] = 0           |
| `test_ti_scaling_sparse`                              | Ti proportional to branch length (2 cases) |
| `test_ti_scaling_dense`                               | Ti proportional to branch length (2 cases) |
| `test_dense_sparse_consistency`                       | Dense ~ sparse for unambiguous sequences   |
| `test_root_state_dense`                               | Root reflects ancestral composition        |
| `test_root_state_sparse`                              | Root reflects Fitch consensus              |
| `test_nij_orientation_multiple_mutations_sparse`      | Multiple mutation directions correct       |
| `test_nij_accumulation_dense`                         | nij accumulates across edges               |
| `test_ti_proportional_to_composition_sparse`          | Ti proportional to state frequency         |
| `test_dense_sparse_nij_direction_agreement`           | Both agree on dominant mutation cell       |
| `test_root_state_total_equals_alignment_length_dense` | sum(root_state) = seq_length               |
| `test_nij_diagonal_zero`                              | nij[i,i] = 0 (dense and sparse)            |
| `test_nij_non_negative`                               | nij[i,j] >= 0 (dense and sparse)           |
| `test_ti_non_negative`                                | Ti[i] >= 0 (dense and sparse)              |

---

## Dense-Sparse Cross-Validation

**File:** [`test_contract_dense_sparse_real.rs`](../../packages/treetime/src/gtr/infer_gtr/__tests__/test_contract_dense_sparse_real.rs)

| Test                                  | Cases | Datasets                                                                        |
| ------------------------------------- | ----- | ------------------------------------------------------------------------------- |
| `test_contract_dense_sparse_real_gtr` | 7     | flu_h3n2_20, ebola_20, rsv_a_20, dengue_20, tb_20, lassa_L_50, mpox_clade_ii_20 |

**Thresholds:**

- pi cosine > 0.997
- W relative Frobenius < 0.21
- mu relative diff < 0.23

---

## Common Functions

**File:** [`test_common.rs`](../../packages/treetime/src/gtr/infer_gtr/__tests__/test_common.rs)

| Test                              | Purpose                                  |
| --------------------------------- | ---------------------------------------- |
| `test_infer_gtr_impl`             | W, pi, mu inference from mutation counts |
| `test_distance_zero_to_uniform`   | Distance: zero -> uniform = 0.5          |
| `test_distance_uniform_to_skewed` | Distance: uniform -> skewed              |
| `test_distance_small_difference`  | Distance: small perturbations            |
| `test_distance_tiny_difference`   | Distance: tiny perturbations             |
| `test_distance_identical`         | Distance: identity = 0                   |

---

## Support Files (no tests)

**File:** [`prop_support.rs`](../../packages/treetime/src/gtr/__tests__/prop_support.rs) - Proptest assertion helpers (`prop_assert_columns_sum_to`, `prop_assert_rows_sum_to`, `prop_assert_detailed_balance`)

**File:** [`test_gtr_numerical_edge_support.rs`](../../packages/treetime/src/gtr/__tests__/test_gtr_numerical_edge/test_gtr_numerical_edge_support.rs) - Stochastic matrix assertion helper (`assert_stochastic_matrix`)
