# Timetree Inference Tests

[Back to index](_index.md)

## Summary

| Category            | Files  | Tests   | Type                 |
| ------------------- | ------ | ------- | -------------------- |
| Coalescent model    | 8      | 52      | Unit + Golden-master |
| Belief propagation  | 4 + 1  | 17      | Unit + Golden-master |
| Convergence         | 2      | 5       | Unit + Integration   |
| Output (confidence) | 1      | 11      | Unit                 |
| Clock filter        | 1      | 4       | Unit                 |
| Relaxed clock       | 1      | 11      | Unit                 |
| Polytomy resolution | 1      | 8       | Unit                 |
| Rerooting           | 1      | 4       | Integration          |
| **Total**           | **20** | **110** | Mixed                |

Test counts include parameterized case expansions (each `#[case]` is one execution).

---

## Coalescent Model Tests

### Piecewise Constant Function

**File:** [`test_piecewise_constant_fn.rs`](../../packages/treetime/src/commands/timetree/coalescent/__tests__/test_piecewise_constant_fn.rs)

| Test                                | Purpose                                    |
| ----------------------------------- | ------------------------------------------ |
| `test_piecewise_constant_eval`      | Step function evaluation at various points |
| `test_piecewise_constant_eval_many` | Batch evaluation matches single-point      |

---

### Piecewise Linear Function

**File:** [`test_piecewise_linear_fn.rs`](../../packages/treetime/src/commands/timetree/coalescent/__tests__/test_piecewise_linear_fn.rs)

| Test                                           | Purpose                               |
| ---------------------------------------------- | ------------------------------------- |
| `test_piecewise_linear_eval_interpolation`     | Linear interpolation within domain    |
| `test_piecewise_linear_eval_extrapolation`     | Constant extrapolation outside domain |
| `test_piecewise_linear_eval_multiple_segments` | Multiple linear segments              |
| `test_piecewise_linear_eval_many`              | Batch evaluation                      |
| `test_piecewise_linear_accessors`              | Getter methods for breakpoints/values |

---

### Event Collection

**File:** [`test_events.rs`](../../packages/treetime/src/commands/timetree/coalescent/__tests__/test_events.rs)

| Test                                           | Purpose                      |
| ---------------------------------------------- | ---------------------------- |
| `test_collect_tree_events_simple`              | Simple tree with 3 children  |
| `test_collect_tree_events_complex`             | Tree with internal nodes     |
| `test_collect_tree_events_sorted`              | Events are time-sorted       |
| `test_collect_tree_events_simultaneous_events` | Multiple events at same time |

**Algorithm:** Event collection from tree topology and dates

---

### Lineage Dynamics

**File:** [`test_lineage_dynamics.rs`](../../packages/treetime/src/commands/timetree/coalescent/__tests__/test_lineage_dynamics.rs)

| Test                                            | Purpose                               |
| ----------------------------------------------- | ------------------------------------- |
| `test_piecewise_constant_eval`                  | Step function evaluation (duplicate)  |
| `test_piecewise_constant_eval_many`             | Batch evaluation (duplicate)          |
| `test_lineage_count_simple_tree`                | Lineage count from simple tree events |
| `test_lineage_count_single_event`               | Single event                          |
| `test_lineage_count_empty_events`               | Error handling for empty input        |
| `test_lineage_count_aggregation`                | Multiple events at same time          |
| `test_lineage_count_decreasing_then_increasing` | Complex event sequence                |
| `test_lineage_count_breakpoints`                | Breakpoint extraction                 |
| `test_lineage_count_negative_deltas`            | Merger events                         |

**Algorithm:** Lineage count dynamics from tree events. File also contains 2 piecewise constant tests.

---

### Integration (Merger Rates)

**File:** [`test_integration.rs`](../../packages/treetime/src/commands/timetree/coalescent/__tests__/test_integration.rs)

| Test                                                         | Purpose                                     |
| ------------------------------------------------------------ | ------------------------------------------- |
| `test_merger_rates`                                          | Basic merger rate computation               |
| `test_merger_rates_edge_cases`                               | Non-integer lineage counts                  |
| `test_merger_rates_large_k`                                  | Large lineage counts                        |
| `test_compute_integral_merger_rate_constant_tc`              | Integration with constant Tc                |
| `test_compute_integral_merger_rate_multiple_segments`        | Multiple segments                           |
| `test_compute_integral_merger_rate_insufficient_points`      | Error handling                              |
| `test_compute_integral_merger_rate_varying_tc`               | Varying Tc (linear)                         |
| `test_compute_integral_merger_rate_varying_tc_many_segments` | Numerical accuracy with analytical solution |

---

### Tc Optimization

**File:** [`test_optimize_tc.rs`](../../packages/treetime/src/commands/timetree/coalescent/__tests__/test_optimize_tc.rs)

| Test                                                 | Cases | Purpose                                               |
| ---------------------------------------------------- | ----- | ----------------------------------------------------- |
| `test_optimize_tc_converges`                         | 1     | Basic convergence                                     |
| `test_optimize_tc_convergence_from_different_starts` | 3     | Convergence from (0.1, 1.0), (0.1, 10.0), (1.0, 10.0) |

4 test executions total. `#[rstest]` with 3 `#[case]` annotations.

---

### Skyline Optimization

**File:** [`test_skyline.rs`](../../packages/treetime/src/commands/timetree/coalescent/__tests__/test_skyline.rs)

| Test                                               | Purpose                              |
| -------------------------------------------------- | ------------------------------------ |
| `test_optimize_skyline_returns_result`             | Basic functionality                  |
| `test_optimize_skyline_tc_distribution_evaluates`  | Tc distribution evaluation           |
| `test_optimize_skyline_log_tc_in_reasonable_range` | Regularization keeps log(Tc) bounded |
| `test_optimize_skyline_larger_tree`                | 8-leaf tree with 10 grid points      |

**Algorithm:** Skyline (piecewise constant Tc) optimization

---

### Golden-Master Coalescent

**File:** [`test_gm_coalescent.rs`](../../packages/treetime/src/commands/timetree/coalescent/__tests__/test_gm_coalescent.rs)

| Test                 | Cases | Datasets and Tc values                                                                                                                                                               |
| -------------------- | ----- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `test_gm_coalescent` | 16    | `flu_h3n2_20` (tc 0.01, 0.1, 1.0, 10.0), `ebola_20` (tc 0.1, 1.0, 10.0), `dengue_20` (tc 1.0, 10.0, 100.0), `rsv_a_20` (tc 1.0, 10.0, 100.0), `mpox_clade_ii_20` (tc 0.1, 1.0, 10.0) |

16 test executions. `#[rstest]` with 16 `#[case]` annotations. Tolerance: 1e-5 max absolute error. Fixtures: `__fixtures__/gm_coalescent_*.json`.

---

## Belief Propagation / Inference Tests

### Backward Pass

**File:** [`test_backward_pass.rs`](../../packages/treetime/src/commands/timetree/inference/__tests__/test_backward_pass.rs)

| Test                                                                  | Purpose                                          |
| --------------------------------------------------------------------- | ------------------------------------------------ |
| `test_backward_pass_computes_internal_node_time`                      | Parent time = child time - branch length         |
| `test_backward_pass_multiplies_child_messages`                        | Message multiplication (constraint intersection) |
| `test_backward_pass_preserves_leaf_time_distribution_with_coalescent` | Coalescent does not overwrite leaf dates         |
| `test_backward_pass_sets_edge_messages`                               | `msg_to_parent` storage                          |
| `test_backward_pass_skips_bad_branch_children`                        | `bad_branch` flag handling                       |
| `test_backward_pass_bad_branch_equivalent_to_removal`                 | `bad_branch` matches tree without that leaf      |

Helper module with `set_leaf_time()` and `set_edge_branch_dist()`.

---

### Runner (Input Mode)

**File:** [`test_runner.rs`](../../packages/treetime/src/commands/timetree/inference/__tests__/test_runner.rs)

| Test                                                           | Purpose                                              |
| -------------------------------------------------------------- | ---------------------------------------------------- |
| `test_create_branch_distributions_input_mode_sets_time_length` | `time_length = branch_length / clock_rate`           |
| `test_input_mode_newick_output_uses_time_lengths`              | Newick serialization uses time lengths               |
| `test_input_mode_gamma_scales_time_length`                     | `time_length = branch_length / (clock_rate * gamma)` |
| `test_input_mode_gamma_default_matches_no_gamma`               | Default gamma=1.0 behavior                           |

### Branch-length likelihood grid

**File:** [`test_branch_length_likelihood.rs`](../../packages/treetime/src/commands/timetree/inference/__tests__/test_branch_length_likelihood.rs)

Direct tests of `compute_branch_length_distribution()`. Empty contribution slices isolate the Poisson indel term; non-empty slices are already exercised by the `test_gm_runner_*` path end-to-end.

| Test                                                                  | Purpose                                                                    |
| --------------------------------------------------------------------- | -------------------------------------------------------------------------- |
| `test_branch_length_likelihood_no_indels_flat_distribution`           | `indel_rate == 0` produces flat distribution (prob == 1 at every sample)   |
| `test_branch_length_likelihood_indel_rate_only_matches_poisson_shape` | `k == 0, mu > 0` follows `exp(-mu (t - t_min))` (shape + peak at grid min) |
| `test_branch_length_likelihood_indel_mle_peak`                        | `k > 0, mu > 0` peaks at Poisson MLE `t_mle = k/mu`                        |
| `test_branch_length_likelihood_indel_mle_peak_with_gamma`             | Gamma compresses the time-domain peak to `t_mle / (clock_rate * gamma)`    |
| `test_branch_length_likelihood_zero_indels_matches_substitution_only` | Explicit no-op check when `indel_count == 0` and `indel_rate == 0`         |
| `test_branch_length_likelihood_rejects_nonpositive_clock_rate`        | Error path: negative clock rate is rejected before the grid is constructed |

---

### Golden-Master Runner Tests

**Directory:** [`test_gm_runner/`](../../packages/treetime/src/commands/timetree/inference/__tests__/test_gm_runner/)

#### Support

**File:** `test_gm_runner_support.rs` - shared fixtures and helpers, no tests.

Contains `DatasetOutputs` struct, `OUTPUTS` lazy static, `load_dates_for_dataset()`, `load_alignment_for_dataset()`. Loaded from `__fixtures__/gm_runner_inputs.json` and `__fixtures__/gm_runner_outputs.json`.

#### Poisson Mode

**File:** `test_gm_runner_poisson.rs`

| Test                     | Cases | Datasets (enabled)        | Tolerance |
| ------------------------ | ----- | ------------------------- | --------- |
| `test_gm_runner_poisson` | 2     | `ebola_20`, `flu_h3n2_20` | 3e-1      |

6 datasets disabled: `dengue_20`, `lassa_L_20`, `mpox_clade_ii_20`, `rsv_a_20`, `tb_20` (zero-length branches), `zika_20` (date column mismatch).

#### Dense Marginal

**File:** `test_gm_runner_marginal_dense.rs`

| Test                            | Cases | Datasets (enabled) | Tolerance |
| ------------------------------- | ----- | ------------------ | --------- |
| `test_gm_runner_marginal_dense` | 1     | `flu_h3n2_20`      | 9e-1      |

`ebola_20` disabled (gap character). Same 6 datasets disabled as poisson.

#### Sparse Marginal

**File:** `test_gm_runner_marginal_sparse.rs`

| Test                             | Cases | Datasets (enabled) | Tolerance |
| -------------------------------- | ----- | ------------------ | --------- |
| `test_gm_runner_marginal_sparse` | 1     | `flu_h3n2_20`      | 9e-1      |

Validates against v0 marginal dense golden values. Same datasets disabled as poisson plus `ebola_20`.

#### ML Branch-Length Pre-Optimization

**File:** `test_gm_runner_pre_optimize.rs`

| Test                                                 | Cases | Datasets      | Purpose                                                |
| ---------------------------------------------------- | ----- | ------------- | ------------------------------------------------------ |
| `test_gm_runner_pre_optimize_changes_branch_lengths` | 1     | `flu_h3n2_20` | Verifies Brent optimization modifies at least one BL   |
| `test_gm_runner_pre_optimize_pipeline_succeeds`      | 1     | `flu_h3n2_20` | Full pipeline with pre-step produces finite node times |

#### Coalescent Integration

**File:** `test_runner_coalescent.rs`

| Test                               | Cases | Tc values      | Purpose                                        |
| ---------------------------------- | ----- | -------------- | ---------------------------------------------- |
| `test_runner_coalescent_completes` | 3     | 0.1, 1.0, 10.0 | Pipeline completes without panic, finite times |

Uses `flu_h3n2_20` dataset. Helper `build_timetree_setup()` in same file.

---

## Convergence Tests

**File:** [`test_metrics.rs`](../../packages/treetime/src/commands/timetree/convergence/__tests__/test_metrics.rs)

| Test                                                     | Purpose                                        |
| -------------------------------------------------------- | ---------------------------------------------- |
| `test_metrics_optimizer_converges_when_n_diff_zero`      | Immediate convergence when `n_diff=0`          |
| `test_metrics_optimizer_continues_when_n_diff_positive`  | Continues until convergence                    |
| `test_metrics_optimizer_stops_at_max_iterations`         | Max iteration limit                            |
| `test_metrics_optimizer_n_resolved_prevents_convergence` | Polytomy resolution prevents early convergence |

**File:** [`test_pipeline.rs`](../../packages/treetime/src/commands/timetree/convergence/__tests__/test_pipeline.rs)

| Test                                 | Purpose                                                                       | Type        |
| ------------------------------------ | ----------------------------------------------------------------------------- | ----------- |
| `test_pipeline_timetree_convergence` | Full pipeline on `flu/h3n2/20`: tracelog CSV, output files, likelihood values | Integration |

---

## Output Tests

**File:** [`test_confidence.rs`](../../packages/treetime/src/commands/timetree/output/__tests__/test_confidence.rs)

| Test                                                                 | Purpose                                  |
| -------------------------------------------------------------------- | ---------------------------------------- |
| `test_extract_confidence_intervals_includes_unnamed_nodes`           | Unnamed nodes included with empty name   |
| `test_extract_confidence_intervals_skips_nodes_without_time`         | Nodes without time excluded              |
| `test_extract_confidence_intervals_uses_date_as_fallback`            | Date used when no distribution           |
| `test_extract_confidence_intervals_with_distribution`                | 90% HPD from distribution                |
| `test_extract_confidence_intervals_sorted_by_key`                    | Sorted by GraphNodeKey (insertion order) |
| `test_combine_confidence_no_contributions`                           | Baseline only                            |
| `test_combine_confidence_single_contribution`                        | One contribution                         |
| `test_combine_confidence_quadrature`                                 | Quadrature sum of deviations             |
| `test_extract_confidence_intervals_clamps_when_date_outside_rate_ci` | Postcondition clamp: date above rate CI  |
| `test_extract_confidence_intervals_clamps_when_date_below_rate_ci`   | Postcondition clamp: date below rate CI  |
| `test_combine_confidence_clipped_to_limits`                          | Clipping to baseline limits              |

**File:** [`test_auspice.rs`](../../packages/treetime/src/commands/timetree/output/__tests__/test_auspice.rs)

| Test                                             | Purpose                                   |
| ------------------------------------------------ | ----------------------------------------- |
| `test_auspice_metadata_structure`                | v2 metadata: panels, colorings, defaults  |
| `test_auspice_root_node_attributes`              | Root node: div=0, num_date, bad_branch    |
| `test_auspice_divergence_from_node_payload`      | div read from NodeTimetree.div            |
| `test_auspice_bad_branch_attribute`              | Yes/No encoding of bad_branch flag        |
| `test_auspice_confidence_intervals_by_key`       | CI lookup by GraphNodeKey                 |
| `test_auspice_unnamed_node_gets_ci_by_key`       | Unnamed node gets CI via key-based lookup |
| `test_auspice_node_without_time_has_no_num_date` | No num_date when time is None             |
| `test_auspice_output_file_is_valid_json`         | File exists and parses as AuspiceTree     |

---

## Optimization Tests

### Clock Filter

**File:** [`test_clock_filter.rs`](../../packages/treetime/src/commands/timetree/optimization/__tests__/test_clock_filter.rs)

| Test                                       | Purpose                                |
| ------------------------------------------ | -------------------------------------- |
| `test_clock_filter_no_outliers_clean_data` | Well-fitting data produces no outliers |
| `test_clock_filter_detects_outlier`        | Extreme deviation detected             |
| `test_clock_filter_iqd_calculation`        | IQD computation                        |
| `test_clock_filter_respects_threshold`     | Higher threshold allows more deviation |

**Algorithm:** IQD-based outlier detection

---

### Relaxed Clock

**File:** [`test_relaxed_clock.rs`](../../packages/treetime/src/commands/timetree/optimization/__tests__/test_relaxed_clock.rs)

| Test                                                         | Purpose                                     |
| ------------------------------------------------------------ | ------------------------------------------- |
| `test_relaxed_clock_default_params_produce_reasonable_gamma` | Gamma in reasonable range                   |
| `test_relaxed_clock_all_gamma_above_minimum`                 | Minimum bound (0.1) enforced                |
| `test_relaxed_clock_uniform_branches_produce_similar_gamma`  | Uniform branches give similar gamma         |
| `test_relaxed_clock_empty_params_uses_defaults`              | Empty params handled gracefully             |
| `test_relaxed_clock_gamma_stored_in_edges`                   | Gamma values stored and differ from default |
| `test_relaxed_clock_high_slack_pulls_toward_one`             | Slack parameter effect                      |
| `test_relaxed_clock_high_coupling_reduces_variation`         | Coupling parameter effect                   |
| `test_relaxed_clock_one_mutation_affects_gamma`              | Different `one_mutation` values matter      |
| `test_relaxed_clock_handles_tiny_one_mutation`               | Defense against near-zero `one_mutation`    |
| `test_relaxed_clock_root_has_branch_penalty`                 | Root uses branch penalty                    |
| `test_relaxed_clock_childless_root_gamma_equals_one`         | Analytical: gamma=1.0 for childless root    |

**Algorithm:** Per-branch rate variation

---

### Polytomy Resolution

**File:** [`test_polytomy.rs`](../../packages/treetime/src/commands/timetree/optimization/__tests__/test_polytomy.rs)

| Test                                                           | Purpose                               |
| -------------------------------------------------------------- | ------------------------------------- |
| `test_find_polytomy_nodes_detects_multifurcation`              | Detection of nodes with >2 children   |
| `test_find_polytomy_nodes_returns_empty_for_binary_tree`       | No false positives on binary tree     |
| `test_resolve_polytomies_reduces_children_count`               | 3-way to 2-way                        |
| `test_resolve_polytomies_no_change_for_binary_tree`            | No changes for binary tree            |
| `test_resolve_polytomies_respects_threshold`                   | High threshold prevents merges        |
| `test_resolve_polytomies_new_node_has_correct_time`            | New node time between parent/children |
| `test_resolve_polytomies_large_polytomy`                       | 5-way polytomy requires 3 merges      |
| `test_prepare_tree_after_topology_change_preserves_leaf_state` | Leaf constraints preserved            |

**Algorithm:** Converting n-way splits to binary

---

### Rerooting

**File:** [`test_reroot.rs`](../../packages/treetime/src/commands/timetree/optimization/__tests__/test_reroot.rs)

| Test                                                 | Purpose                         | Type        |
| ---------------------------------------------------- | ------------------------------- | ----------- |
| `test_reroot_tree_sparse_with_edge_split`            | Reroot with sparse partition    | Integration |
| `test_sparse_reroot_inverts_subs_and_indels_on_path` | Substitution/indel inversion    | Unit        |
| `test_sparse_reroot_inverts_edge_mutations`          | Detailed mutation inversion     | Unit        |
| `test_reroot_tree_sparse_flow_does_not_panic`        | Two successive reroots complete | Regression  |

---

## Module Structure

Each `__tests__/mod.rs` contains only `mod` declarations.

| Directory                                   | Module declarations                                                                                                                                                            |
| ------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `coalescent/__tests__/mod.rs`               | `test_events`, `test_gm_coalescent`, `test_integration`, `test_lineage_dynamics`, `test_optimize_tc`, `test_piecewise_constant_fn`, `test_piecewise_linear_fn`, `test_skyline` |
| `inference/__tests__/mod.rs`                | `test_backward_pass`, `test_gm_runner`, `test_runner`                                                                                                                          |
| `inference/__tests__/test_gm_runner/mod.rs` | `test_gm_runner_marginal_dense`, `test_gm_runner_marginal_sparse`, `test_gm_runner_poisson`, `test_gm_runner_pre_optimize`, `test_gm_runner_support`, `test_runner_coalescent` |
| `convergence/__tests__/mod.rs`              | `test_metrics`, `test_pipeline`                                                                                                                                                |
| `output/__tests__/mod.rs`                   | `test_confidence`                                                                                                                                                              |
| `optimization/__tests__/mod.rs`             | `test_clock_filter`, `test_polytomy`, `test_relaxed_clock`, `test_reroot`                                                                                                      |

---

## Known Test Limitations

1. **Golden-master datasets disabled**: Multiple datasets disabled due to zero-length branches:
   `dengue_20`, `lassa_L_20`, `mpox_clade_ii_20`, `rsv_a_20`, `tb_20`, `zika_20`

2. **Tolerance variations**: Root node dominates max diff in Poisson mode tests. Non-root nodes typically agree within 1e-2.

3. **Dense marginal**: `ebola_20` disabled due to gap character handling in alphabet.

4. **Duplicate tests**: `test_lineage_dynamics.rs` contains 2 piecewise constant tests that duplicate those in `test_piecewise_constant_fn.rs`.
