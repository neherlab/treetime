# Mugration Test Coverage

[Back to index](_index.md)

## Summary

| Category                                                        | Type                                           |
| --------------------------------------------------------------- | ---------------------------------------------- |
| [Golden master (v0 parity)](#golden-master-v0-parity)           | Golden-master                                  |
| [Structural / unit](#structural--unit)                          | Unit                                           |
| [Algorithm invariants](#algorithm-invariants)                   | Unit                                           |
| [Discrete marginal](#discrete-marginal)                         | Unit                                           |
| [Comment output](#comment-output)                               | Unit                                           |
| [Brent optimizer](#brent-optimizer)                             | Unit                                           |
| [Partition / discrete](representation.md#normalize-discrete-1d) | Unit (see [Representation](representation.md)) |

---

## Golden Master (v0 parity)

**Test:** [`packages/treetime/src/commands/mugration/__tests__/test_gm_mugration.rs`](../../packages/treetime/src/commands/mugration/__tests__/test_gm_mugration.rs)

**Impl:**

- [`packages/treetime/src/commands/mugration/run.rs`](../../packages/treetime/src/commands/mugration/run.rs)
- [`packages/treetime/src/commands/mugration/input.rs`](../../packages/treetime/src/commands/mugration/input.rs)
- [`packages/treetime/src/commands/mugration/output.rs`](../../packages/treetime/src/commands/mugration/output.rs)

| Test                                      | Datasets                                   | Notes                                                                            |
| ----------------------------------------- | ------------------------------------------ | -------------------------------------------------------------------------------- |
| `test_gm_mugration_outputs`               | zika_20_country, zika_20_country_weights   | passing                                                                          |
| `test_gm_mugration_outputs_v1_divergence` | lassa, dengue, tb, rsv, mpox               | **ignored**: intentional v1 improvements (pseudo-count pi, root-state filtering) |
| `test_gm_mugration_confidence_zika`       | zika                                       | passing, **2e-2** tolerance                                                      |
| `test_gm_mugration_confidence_outputs`    | zika_weights, lassa, dengue, tb, rsv, mpox | **ignored**: confidence profiles exceed 1e-2 at multiple nodes                   |

---

## Structural / Unit

**Test:** [`packages/treetime/src/commands/mugration/__tests__/test_run.rs`](../../packages/treetime/src/commands/mugration/__tests__/test_run.rs)

**Impl:** [`packages/treetime/src/commands/mugration/run.rs`](../../packages/treetime/src/commands/mugration/run.rs)

| Test                                                                    | Purpose                                                       |
| ----------------------------------------------------------------------- | ------------------------------------------------------------- |
| `test_run_validate_weight_coverage_rejects_above_threshold`             | Weight coverage validation rejects above threshold            |
| `test_run_validate_weight_coverage_accepts_at_threshold`                | Weight coverage validation accepts at threshold               |
| `test_run_validate_weight_coverage_excludes_missing_data_marker`        | Missing data marker excluded from coverage                    |
| `test_run_validate_weight_coverage_full_coverage`                       | Full coverage passes validation                               |
| `test_run_compute_pi_from_weights_normalizes`                           | Pi from weights is normalized                                 |
| `test_run_compute_pi_from_weights_uses_mean_for_missing`                | Missing weights use mean value                                |
| `test_run_compute_pi_uniform`                                           | Uniform pi computation                                        |
| `test_run_apply_pseudo_counts_with_value`                               | Pseudo-count effect on pi                                     |
| `test_run_apply_pseudo_counts_without_value`                            | No pseudo-count preserves pi                                  |
| `test_run_apply_pseudo_counts_preserves_normalization`                  | Pseudo-counts preserve normalization                          |
| `test_execute_mugration_simple_tree`                                    | Basic 2-state tree: states, assignments, confidence structure |
| `test_execute_mugration_with_weights`                                   | Weight-based pi computation, weighted reconstruction          |
| `test_execute_mugration_with_weights_includes_unobserved_weight_states` | Unobserved weight states included in alphabet                 |
| `test_execute_mugration_with_pseudo_counts`                             | Pseudo-count effect on pi with weights                        |
| `test_execute_mugration_sampling_bias_correction`                       | mu scaled by correction factor                                |
| `test_execute_mugration_rejects_single_state`                           | Error on < 2 states                                           |

---

## Algorithm Invariants

**Test:** [`packages/treetime/src/commands/mugration/__tests__/test_run.rs`](../../packages/treetime/src/commands/mugration/__tests__/test_run.rs) (same file)

**Impl:** [`packages/treetime/src/commands/mugration/run.rs`](../../packages/treetime/src/commands/mugration/run.rs)

| Test                                           | Purpose                                            |
| ---------------------------------------------- | -------------------------------------------------- |
| `test_iterative_refinement_changes_model`      | mu/pi differ between iterations=0 and iterations=5 |
| `test_iterative_refinement_pi_reflects_data`   | pi shifts toward observed frequencies with free pi |
| `test_zero_iterations_preserves_initial_model` | Basic parameter validity at iterations=0           |

---

## Discrete Marginal

**Test:** [`packages/treetime/src/commands/mugration/__tests__/test_discrete_marginal.rs`](../../packages/treetime/src/commands/mugration/__tests__/test_discrete_marginal.rs)

**Impl:** [`packages/treetime/src/commands/mugration/discrete_marginal.rs`](../../packages/treetime/src/commands/mugration/discrete_marginal.rs)

| Test                                                                               | Purpose                                          |
| ---------------------------------------------------------------------------------- | ------------------------------------------------ |
| `test_discrete_marginal_attach_traits_maps_observed_and_missing_profiles`          | Observed states get one-hot, missing get uniform |
| `test_discrete_marginal_attach_traits_rejects_tree_leaf_missing_from_metadata`     | Tree leaf missing from metadata rejected         |
| `test_discrete_marginal_attach_traits_rejects_metadata_name_missing_from_tree`     | Metadata name missing from tree rejected         |
| `test_discrete_marginal_passes_normalize_backward_and_forward_profiles`            | Backward and forward profiles normalized         |
| `test_discrete_marginal_run_returns_finite_log_lh_and_reconstructs_internal_trait` | Finite log-LH and internal trait reconstruction  |

---

## Comment Output

**Test:** [`packages/treetime/src/commands/mugration/__tests__/test_comment_output.rs`](../../packages/treetime/src/commands/mugration/__tests__/test_comment_output.rs)

**Impl:** [`packages/treetime/src/commands/mugration/comment_provider.rs`](../../packages/treetime/src/commands/mugration/comment_provider.rs)

| Test                                               | Purpose                                         |
| -------------------------------------------------- | ----------------------------------------------- |
| `test_mugration_annotated_tree_has_trait_comments` | Nexus annotation format includes trait comments |

---

## Brent Optimizer

**Test:** [`packages/treetime/src/commands/mugration/gtr_refinement.rs`](../../packages/treetime/src/commands/mugration/gtr_refinement.rs) (inline `#[cfg(test)]`)

**Impl:** [`packages/treetime/src/commands/mugration/gtr_refinement.rs`](../../packages/treetime/src/commands/mugration/gtr_refinement.rs)

| Test                                                 | Purpose                                        |
| ---------------------------------------------------- | ---------------------------------------------- |
| `test_brent_minimize_bracketed_quadratic`            | Minimum of (x-3)^2 with bracket triple         |
| `test_brent_minimize_bracketed_shifted_quadratic`    | Minimum of (x-0.7)^2+1 with bracket triple     |
| `test_brent_minimize_bracketed_cosine`               | Minimum of cos(x) on [2,5] with bracket triple |
| `test_brent_minimize_bracketed_starts_near_interior` | Fast convergence from near-optimal interior    |

---

## Partition / Discrete

Moved to [Representation Tests: Normalize (Discrete 1D)](representation.md#normalize-discrete-1d).
