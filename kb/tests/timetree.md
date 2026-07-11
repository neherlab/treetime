# Timetree Inference Tests

[Back to index](README.md)

## Summary

| Category                                                   | Type                 |
| ---------------------------------------------------------- | -------------------- |
| [Coalescent model](#coalescent-model-tests)                | Unit + Golden-master |
| [Total log-likelihood](#total-log-likelihood)              | Unit + Golden-master |
| [Time coordinate](#time-coordinate)                        | Unit                 |
| [Belief propagation](#belief-propagation--inference-tests) | Unit + Golden-master |
| [Convergence](#convergence-tests)                          | Unit + Integration   |
| [Output (confidence)](#output-tests)                       | Unit                 |
| [Output (Auspice)](#auspice-output)                        | Unit                 |
| [Clock filter](#clock-filter)                              | Unit                 |
| [Relaxed clock](#relaxed-clock)                            | Unit                 |
| [Polytomy resolution](#polytomy-resolution)                | Unit                 |
| [Rerooting](#rerooting)                                    | Integration          |

Parameterized tests use `#[case]` expansions via rstest.

---

## Coalescent Model Tests

### Piecewise Constant Function

**Test:** [`packages/treetime-grid/src/__tests__/piecewise_constant_fn.rs`](../../packages/treetime-grid/src/__tests__/piecewise_constant_fn.rs)

**Impl:** [`packages/treetime-grid/src/piecewise_constant_fn.rs`](../../packages/treetime-grid/src/piecewise_constant_fn.rs)

| Test                                | Purpose                                    |
| ----------------------------------- | ------------------------------------------ |
| `test_piecewise_constant_eval`      | Step function evaluation at various points |
| `test_piecewise_constant_eval_many` | Batch evaluation matches single-point      |

---

### Piecewise Linear Function

**Test:** [`packages/treetime-grid/src/__tests__/piecewise_linear_fn.rs`](../../packages/treetime-grid/src/__tests__/piecewise_linear_fn.rs)

**Impl:** [`packages/treetime-grid/src/piecewise_linear_fn.rs`](../../packages/treetime-grid/src/piecewise_linear_fn.rs)

| Test                                           | Purpose                               |
| ---------------------------------------------- | ------------------------------------- |
| `test_piecewise_linear_eval_interpolation`     | Linear interpolation within domain    |
| `test_piecewise_linear_eval_extrapolation`     | Constant extrapolation outside domain |
| `test_piecewise_linear_eval_multiple_segments` | Multiple linear segments              |
| `test_piecewise_linear_eval_many`              | Batch evaluation                      |
| `test_piecewise_linear_accessors`              | Getter methods for breakpoints/values |

---

### Event Collection

**Test:** [`packages/treetime/src/coalescent/__tests__/test_events.rs`](../../packages/treetime/src/coalescent/__tests__/test_events.rs)

**Impl:** [`packages/treetime/src/coalescent/events.rs`](../../packages/treetime/src/coalescent/events.rs)

| Test                                           | Purpose                      |
| ---------------------------------------------- | ---------------------------- |
| `test_collect_tree_events_simple`              | Simple tree with 3 children  |
| `test_collect_tree_events_complex`             | Tree with internal nodes     |
| `test_collect_tree_events_sorted`              | Events are time-sorted       |
| `test_collect_tree_events_simultaneous_events` | Multiple events at same time |

---

### Lineage Dynamics

**Test:** [`packages/treetime/src/coalescent/__tests__/test_lineage_dynamics.rs`](../../packages/treetime/src/coalescent/__tests__/test_lineage_dynamics.rs)

**Impl:** [`packages/treetime/src/coalescent/lineage_dynamics.rs`](../../packages/treetime/src/coalescent/lineage_dynamics.rs)

| Test                                            | Purpose                                                       |
| ----------------------------------------------- | ------------------------------------------------------------- |
| `test_piecewise_constant_eval`                  | Step function evaluation (duplicate of piecewise_constant_fn) |
| `test_piecewise_constant_eval_many`             | Batch evaluation (duplicate of piecewise_constant_fn)         |
| `test_lineage_count_simple_tree`                | Lineage count from simple tree events                         |
| `test_lineage_count_single_event`               | Single event                                                  |
| `test_lineage_count_empty_events`               | Error handling for empty input                                |
| `test_lineage_count_aggregation`                | Multiple events at same time                                  |
| `test_lineage_count_decreasing_then_increasing` | Complex event sequence                                        |
| `test_lineage_count_breakpoints`                | Breakpoint extraction                                         |
| `test_lineage_count_negative_deltas`            | Merger events                                                 |

---

### Integration (Merger Rates)

**Test:** [`packages/treetime/src/coalescent/__tests__/test_integration.rs`](../../packages/treetime/src/coalescent/__tests__/test_integration.rs)

**Impl:** [`packages/treetime/src/coalescent/integration.rs`](../../packages/treetime/src/coalescent/integration.rs)

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

**Test:** [`packages/treetime/src/coalescent/__tests__/test_optimize_tc.rs`](../../packages/treetime/src/coalescent/__tests__/test_optimize_tc.rs)

**Impl:** [`packages/treetime/src/coalescent/optimize_tc.rs`](../../packages/treetime/src/coalescent/optimize_tc.rs)

| Test                                                 | Purpose                                               |
| ---------------------------------------------------- | ----------------------------------------------------- |
| `test_optimize_tc_converges`                         | Basic convergence                                     |
| `test_optimize_tc_convergence_from_different_starts` | Convergence from (0.1, 1.0), (0.1, 10.0), (1.0, 10.0) |

---

### Skyline Optimization

**Test:** [`packages/treetime/src/coalescent/__tests__/test_skyline.rs`](../../packages/treetime/src/coalescent/__tests__/test_skyline.rs)

**Impl:** [`packages/treetime/src/coalescent/skyline.rs`](../../packages/treetime/src/coalescent/skyline.rs)

| Test                                               | Purpose                              |
| -------------------------------------------------- | ------------------------------------ |
| `test_optimize_skyline_returns_result`             | Basic functionality                  |
| `test_optimize_skyline_tc_distribution_evaluates`  | Tc distribution evaluation           |
| `test_optimize_skyline_log_tc_in_reasonable_range` | Regularization keeps log(Tc) bounded |
| `test_optimize_skyline_larger_tree`                | 8-leaf tree with 10 grid points      |

---

### Golden-Master Coalescent

**Test:** [`packages/treetime/src/coalescent/__tests__/test_gm_coalescent.rs`](../../packages/treetime/src/coalescent/__tests__/test_gm_coalescent.rs)

**Impl:** [`packages/treetime/src/coalescent/coalescent.rs`](../../packages/treetime/src/coalescent/coalescent.rs)

| Test                 | Datasets and Tc values                                                                                                                                                               |
| -------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `test_gm_coalescent` | `flu_h3n2_20` (tc 0.01, 0.1, 1.0, 10.0), `ebola_20` (tc 0.1, 1.0, 10.0), `dengue_20` (tc 1.0, 10.0, 100.0), `rsv_a_20` (tc 1.0, 10.0, 100.0), `mpox_clade_ii_20` (tc 0.1, 1.0, 10.0) |

Tolerance: 1e-5 max absolute error. Fixtures: [`gm_coalescent_*.json`](../../packages/treetime/src/coalescent/__tests__/__fixtures__/).

---

## Total Log-Likelihood

**Test:** [`packages/treetime/src/coalescent/__tests__/test_total_lh.rs`](../../packages/treetime/src/coalescent/__tests__/test_total_lh.rs)

**Impl:** [`packages/treetime/src/coalescent/total_lh.rs`](../../packages/treetime/src/coalescent/total_lh.rs)

| Test                                           | Purpose                                                      |
| ---------------------------------------------- | ------------------------------------------------------------ |
| `test_total_lh_returns_finite_value`           | Coalescent total LH is finite                                |
| `test_total_lh_negative_for_reasonable_tc`     | Coalescent log-likelihood is negative for non-trivial trees  |
| `test_total_lh_finite_for_tc` (4 cases)        | LH finite across range of Tc values (0.1, 1, 10, 100)        |
| `test_total_lh_differs_across_tc_values`       | Different Tc values produce different LH                     |
| `test_total_lh_monotonic_near_optimum`         | LH peaks near optimal Tc                                     |
| `test_total_lh_matches_optimize_tc_likelihood` | compute_coalescent_total_lh matches optimize_tc output       |
| `test_total_lh_with_formula_distribution`      | Constant Formula distribution matches Distribution::constant |

**Test:** [`packages/treetime/src/coalescent/__tests__/test_gm_total_lh.rs`](../../packages/treetime/src/coalescent/__tests__/test_gm_total_lh.rs)

**Impl:** [`packages/treetime/src/coalescent/total_lh.rs`](../../packages/treetime/src/coalescent/total_lh.rs)

| Test                        | Purpose                                                            |
| --------------------------- | ------------------------------------------------------------------ |
| `test_gm_total_lh_binary`   | Golden master: binary tree total LH matches v0 at four Tc values   |
| `test_gm_total_lh_polytomy` | Golden master: polytomy tree total LH matches v0 at four Tc values |

---

## Time Coordinate

**Test:** [`packages/treetime/src/coalescent/time_coordinate.rs`](../../packages/treetime/src/coalescent/time_coordinate.rs) (inline `#[cfg(test)]`)

**Impl:** [`packages/treetime/src/coalescent/time_coordinate.rs`](../../packages/treetime/src/coalescent/time_coordinate.rs)

| Test                             | Purpose                                                       |
| -------------------------------- | ------------------------------------------------------------- |
| `test_calendar_to_tbp_roundtrip` | CalendarTime -> Tbp -> CalendarTime roundtrip preserves value |
| `test_tbp_at_present_is_zero`    | Present time converts to Tbp(0)                               |
| `test_tbp_arithmetic`            | Tbp subtraction and addition produce correct values           |
| `test_calendar_time_max`         | CalendarTime::max returns the later time                      |

---

## Belief Propagation / Inference Tests

### Backward Pass

**Test:** [`packages/treetime/src/timetree/inference/__tests__/test_backward_pass.rs`](../../packages/treetime/src/timetree/inference/__tests__/test_backward_pass.rs)

**Impl:** [`packages/treetime/src/timetree/inference/backward_pass.rs`](../../packages/treetime/src/timetree/inference/backward_pass.rs)

| Test                                                                    | Purpose                                          |
| ----------------------------------------------------------------------- | ------------------------------------------------ |
| `test_backward_pass_computes_internal_node_time`                        | Parent time = child time - branch length         |
| `test_backward_pass_multiplies_child_messages`                          | Message multiplication (constraint intersection) |
| `test_backward_pass_preserves_leaf_time_distribution_with_coalescent`   | Coalescent does not overwrite leaf dates         |
| `test_backward_pass_preserves_internal_time_with_large_coalescent_cost` | Large coalescent costs do not underflow to Empty  |
| `test_backward_pass_sets_edge_messages`                                 | `msg_to_parent` storage                           |
| `test_backward_pass_skips_bad_branch_children`                          | `bad_branch` flag handling                       |
| `test_backward_pass_bad_branch_equivalent_to_removal`                   | `bad_branch` matches tree without that leaf       |

---

### Runner (Input Mode)

**Test:** [`packages/treetime/src/timetree/inference/__tests__/test_runner.rs`](../../packages/treetime/src/timetree/inference/__tests__/test_runner.rs)

**Impl:** [`packages/treetime/src/timetree/inference/runner.rs`](../../packages/treetime/src/timetree/inference/runner.rs)

| Test                                                           | Purpose                                              |
| -------------------------------------------------------------- | ---------------------------------------------------- |
| `test_create_branch_distributions_input_mode_sets_time_length` | `time_length = branch_length / clock_rate`           |
| `test_input_mode_newick_output_uses_time_lengths`              | Newick serialization uses time lengths               |
| `test_input_mode_gamma_scales_time_length`                     | `time_length = branch_length / (clock_rate * gamma)` |
| `test_input_mode_gamma_default_matches_no_gamma`               | Default gamma=1.0 behavior                           |

### Branch-Length Likelihood Grid

**Test:** [`packages/treetime/src/timetree/inference/__tests__/test_branch_length_likelihood.rs`](../../packages/treetime/src/timetree/inference/__tests__/test_branch_length_likelihood.rs)

**Impl:** [`packages/treetime/src/timetree/inference/branch_length_likelihood.rs`](../../packages/treetime/src/timetree/inference/branch_length_likelihood.rs)

| Test                                                                  | Purpose                                                            | Notes            |
| --------------------------------------------------------------------- | ------------------------------------------------------------------ | ---------------- |
| `test_branch_length_likelihood_no_indels_flat_distribution`           | `indel_rate == 0` produces flat distribution                       |                  |
| `test_branch_length_likelihood_indel_rate_only_matches_poisson_shape` | `k == 0, mu > 0` follows `exp(-mu (t - t_min))`                    | **1e-2** epsilon |
| `test_branch_length_likelihood_indel_mle_peak`                        | `k > 0, mu > 0` peaks at Poisson MLE `t_mle = k/mu`                |                  |
| `test_branch_length_likelihood_indel_mle_peak_with_gamma`             | Gamma compresses the time-domain peak                              |                  |
| `test_branch_length_likelihood_zero_indels_matches_substitution_only` | Explicit no-op check when `indel_count == 0` and `indel_rate == 0` |                  |
| `test_branch_length_likelihood_rejects_nonpositive_clock_rate`        | Error path: negative clock rate rejected                           |                  |

---

### Golden-Master Runner Tests

**Directory:** [`test_gm_runner/`](../../packages/treetime/src/timetree/inference/__tests__/test_gm_runner/)

#### Support

**Test:** [`packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_support.rs`](../../packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_support.rs) - shared fixtures and helpers, no tests.

#### Poisson Mode

**Test:** [`packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_poisson.rs`](../../packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_poisson.rs)

**Impl:**

- [`packages/treetime/src/timetree/inference/backward_pass.rs`](../../packages/treetime/src/timetree/inference/backward_pass.rs)
- [`packages/treetime/src/timetree/inference/forward_pass.rs`](../../packages/treetime/src/timetree/inference/forward_pass.rs)
- [`packages/treetime/src/timetree/inference/runner.rs`](../../packages/treetime/src/timetree/inference/runner.rs)

| Test                     | Datasets              | Tolerance | Notes                                                                      |
| ------------------------ | --------------------- | --------- | -------------------------------------------------------------------------- |
| `test_gm_runner_poisson` | ebola_20, flu_h3n2_20 | **3e-1**  | 6 datasets disabled (see [\README.md](README.md#commented-out-test-cases)) |

#### Dense Marginal

**Test:** [`packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs`](../../packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs)

**Impl:** [`packages/treetime/src/timetree/inference/runner.rs`](../../packages/treetime/src/timetree/inference/runner.rs)

| Test                            | Datasets    | Tolerance | Notes                                                                    |
| ------------------------------- | ----------- | --------- | ------------------------------------------------------------------------ |
| `test_gm_runner_marginal_dense` | flu_h3n2_20 | **9e-1**  | **ignored**: golden master datasets not yet passing. 7 datasets disabled |

#### Sparse Marginal

**Test:** [`packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_sparse.rs`](../../packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_sparse.rs)

**Impl:** [`packages/treetime/src/timetree/inference/runner.rs`](../../packages/treetime/src/timetree/inference/runner.rs)

| Test                             | Datasets    | Tolerance | Notes               |
| -------------------------------- | ----------- | --------- | ------------------- |
| `test_gm_runner_marginal_sparse` | flu_h3n2_20 | **9e-1**  | 6 datasets disabled |

#### ML Branch-Length Pre-Optimization

**Test:** [`packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_pre_optimize.rs`](../../packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_gm_runner_pre_optimize.rs)

**Impl:** [`packages/treetime/src/timetree/inference/runner.rs`](../../packages/treetime/src/timetree/inference/runner.rs)

| Test                                                 | Purpose                                                |
| ---------------------------------------------------- | ------------------------------------------------------ |
| `test_gm_runner_pre_optimize_changes_branch_lengths` | Verifies Brent optimization modifies at least one BL   |
| `test_gm_runner_pre_optimize_pipeline_succeeds`      | Full pipeline with pre-step produces finite node times |

#### Coalescent Integration

**Test:** [`packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_runner_coalescent.rs`](../../packages/treetime/src/timetree/inference/__tests__/test_gm_runner/test_runner_coalescent.rs)

**Impl:** [`packages/treetime/src/timetree/inference/runner.rs`](../../packages/treetime/src/timetree/inference/runner.rs)

| Test                               | Datasets                        | Notes                                               |
| ---------------------------------- | ------------------------------- | --------------------------------------------------- |
| `test_runner_coalescent_completes` | flu_h3n2_20 (Tc 0.1, 1.0, 10.0) | **ignored**: golden master datasets not yet passing |

---

## Convergence Tests

**Test:** [`packages/treetime/src/timetree/convergence/__tests__/test_optimizer.rs`](../../packages/treetime/src/timetree/convergence/__tests__/test_optimizer.rs)

**Impl:** [`packages/treetime/src/timetree/convergence/metrics.rs`](../../packages/treetime/src/timetree/convergence/metrics.rs)

| Test                                                     | Purpose                                        |
| -------------------------------------------------------- | ---------------------------------------------- |
| `test_metrics_optimizer_converges_when_n_diff_zero`      | Immediate convergence when `n_diff=0`          |
| `test_metrics_optimizer_continues_when_n_diff_positive`  | Continues until convergence                    |
| `test_metrics_optimizer_stops_at_max_iterations`         | Max iteration limit                            |
| `test_metrics_optimizer_n_resolved_prevents_convergence` | Polytomy resolution prevents early convergence |

**Test:** [`packages/treetime/src/timetree/convergence/__tests__/test_optimizer.rs`](../../packages/treetime/src/timetree/convergence/__tests__/test_optimizer.rs)

**Impl:** [`packages/treetime/src/commands/timetree/run.rs`](../../packages/treetime/src/commands/timetree/run.rs)

| Test                                 | Purpose                                                                       | Type        |
| ------------------------------------ | ----------------------------------------------------------------------------- | ----------- |
| `test_pipeline_timetree_convergence` | Full pipeline on `flu/h3n2/20`: tracelog CSV, output files, likelihood values | Integration |

---

## Output Tests

### Confidence Extraction

**Test:** [`packages/treetime/src/commands/timetree/output/__tests__/test_confidence_extract.rs`](../../packages/treetime/src/commands/timetree/output/__tests__/test_confidence_extract.rs)

**Impl:** [`packages/treetime/src/timetree/confidence.rs`](../../packages/treetime/src/timetree/confidence.rs)

| Test                                                                 | Purpose                                            | Notes            |
| -------------------------------------------------------------------- | -------------------------------------------------- | ---------------- |
| `test_extract_confidence_intervals_includes_unnamed_nodes`           | Unnamed nodes included with empty name             |                  |
| `test_extract_confidence_intervals_skips_nodes_without_time`         | Nodes without time excluded                        |                  |
| `test_extract_confidence_intervals_uses_date_as_fallback`            | Date used when no distribution                     |                  |
| `test_extract_confidence_intervals_with_distribution`                | 90% HPD from distribution                          |                  |
| `test_extract_confidence_intervals_sorted_by_key`                    | Sorted by GraphNodeKey (insertion order)           |                  |
| `test_extract_confidence_intervals_rate_only`                        | Rate susceptibility dates produce z-score based CI |                  |
| `test_extract_confidence_intervals_combined_wider_than_either`       | Quadrature combination wider than either source    |                  |
| `test_extract_confidence_intervals_clamps_when_date_outside_rate_ci` | Postcondition clamp: date above rate CI            | **1e-3** epsilon |
| `test_extract_confidence_intervals_clamps_when_date_below_rate_ci`   | Postcondition clamp: date below rate CI            | **1e-3** epsilon |
| `test_extract_confidence_intervals_skewed_distribution_hpd`          | HPD region for skewed distribution                 |                  |

### Confidence Combination

**Test:** [`packages/treetime/src/commands/timetree/output/__tests__/test_confidence_combine.rs`](../../packages/treetime/src/commands/timetree/output/__tests__/test_confidence_combine.rs)

**Impl:** [`packages/treetime/src/timetree/confidence.rs`](../../packages/treetime/src/timetree/confidence.rs)

| Test                                          | Purpose                      |
| --------------------------------------------- | ---------------------------- |
| `test_combine_confidence_no_contributions`    | Baseline only                |
| `test_combine_confidence_single_contribution` | One contribution             |
| `test_combine_confidence_quadrature`          | Quadrature sum of deviations |
| `test_combine_confidence_clipped_to_limits`   | Clipping to baseline limits  |

### Rate Uncertainty

**Test:** [`packages/treetime/src/commands/timetree/output/__tests__/test_confidence_rate.rs`](../../packages/treetime/src/commands/timetree/output/__tests__/test_confidence_rate.rs)

**Impl:** [`packages/treetime/src/timetree/confidence.rs`](../../packages/treetime/src/timetree/confidence.rs)

| Test                                                            | Purpose                                             |
| --------------------------------------------------------------- | --------------------------------------------------- |
| `test_quantile_to_zscore` (7 cases)                             | Probit function maps quantiles to z-scores          |
| `test_date_uncertainty_due_to_rate` (5 cases)                   | Date uncertainty from rate susceptibility           |
| `test_determine_rate_std_explicit_clock_std_dev`                | Explicit --clock-std-dev passed through             |
| `test_determine_rate_std_rejects_negative`                      | Negative --clock-std-dev rejected                   |
| `test_determine_rate_std_rejects_zero`                          | Zero --clock-std-dev rejected                       |
| `test_determine_rate_std_none_without_covariation`              | None returned when no covariation                   |
| `test_determine_rate_std_from_covariance_matrix`                | Rate std derived from covariance matrix             |
| `test_determine_rate_std_none_for_fixed_clock_with_covariation` | None returned for fixed clock even with covariation |

### Auspice Output

**Test:** [`packages/treetime/src/commands/timetree/output/__tests__/test_auspice.rs`](../../packages/treetime/src/commands/timetree/output/__tests__/test_auspice.rs)

**Impl:** [`packages/treetime/src/commands/timetree/output/auspice.rs`](../../packages/treetime/src/commands/timetree/output/auspice.rs)

| Test                                             | Purpose                                           |
| ------------------------------------------------ | ------------------------------------------------- |
| `test_auspice_metadata_structure`                | v2 metadata: panels, colorings, defaults          |
| `test_auspice_root_node_attributes`              | Root node: div=0, num_date, bad_branch            |
| `test_auspice_divergence_from_node_payload`      | div read from NodeTimetree.div                    |
| `test_auspice_bad_branch_attribute`              | Yes/No encoding of bad_branch flag                |
| `test_auspice_confidence_intervals_by_key`       | CI lookup by GraphNodeKey                         |
| `test_auspice_unnamed_node_gets_ci_by_key`       | Unnamed node gets CI via key-based lookup         |
| `test_auspice_node_without_time_has_no_num_date` | No num_date when time is None                     |
| `test_auspice_rejects_nan_div`                   | NaN div produces error mentioning node name       |
| `test_auspice_rejects_infinite_time`             | Infinite time produces error mentioning node name |
| `test_auspice_output_file_is_valid_json`         | File exists and parses as AuspiceTree             |

---

## Optimization Tests

### Clock Filter

**Test:** [`packages/treetime/src/timetree/optimization/__tests__/test_clock_filter.rs`](../../packages/treetime/src/timetree/optimization/__tests__/test_clock_filter.rs)

**Impl:** [`packages/treetime/src/clock/clock_filter.rs`](../../packages/treetime/src/clock/clock_filter.rs)

| Test                                       | Purpose                                |
| ------------------------------------------ | -------------------------------------- |
| `test_clock_filter_no_outliers_clean_data` | Well-fitting data produces no outliers |
| `test_clock_filter_detects_outlier`        | Extreme deviation detected             |
| `test_clock_filter_iqd_calculation`        | IQD computation                        |
| `test_clock_filter_respects_threshold`     | Higher threshold allows more deviation |

---

### Relaxed Clock

**Test:** [`packages/treetime/src/timetree/optimization/__tests__/test_relaxed_clock.rs`](../../packages/treetime/src/timetree/optimization/__tests__/test_relaxed_clock.rs)

**Impl:** [`packages/treetime/src/timetree/optimization/relaxed_clock.rs`](../../packages/treetime/src/timetree/optimization/relaxed_clock.rs)

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

---

### Polytomy Resolution

**Test:** [`packages/treetime/src/timetree/optimization/__tests__/test_polytomy.rs`](../../packages/treetime/src/timetree/optimization/__tests__/test_polytomy.rs)

**Impl:** [`packages/treetime/src/timetree/optimization/polytomy.rs`](../../packages/treetime/src/timetree/optimization/polytomy.rs)

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

---

### Rerooting

**Test:** [`packages/treetime/src/timetree/optimization/__tests__/test_reroot.rs`](../../packages/treetime/src/timetree/optimization/__tests__/test_reroot.rs)

**Impl:** [`packages/treetime/src/timetree/optimization/reroot.rs`](../../packages/treetime/src/timetree/optimization/reroot.rs)

| Test                                                 | Purpose                         | Type        |
| ---------------------------------------------------- | ------------------------------- | ----------- |
| `test_reroot_tree_sparse_with_edge_split`            | Reroot with sparse partition    | Integration |
| `test_sparse_reroot_inverts_subs_and_indels_on_path` | Substitution/indel inversion    | Unit        |
| `test_sparse_reroot_inverts_edge_mutations`          | Detailed mutation inversion     | Unit        |
| `test_reroot_tree_sparse_flow_does_not_panic`        | Two successive reroots complete | Regression  |

---

## Module Structure

Each `__tests__/mod.rs` contains only `mod` declarations.

| Directory                                                                                                                              | Module declarations                                                                                                                                                                                                            |
| -------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| [`coalescent/__tests__/mod.rs`](../../packages/treetime/src/coalescent/__tests__/mod.rs)                             | `helpers`, `test_events`, `test_gm_coalescent`, `test_gm_total_lh`, `test_integration`, `test_lineage_dynamics`, `test_optimize_tc`, `test_piecewise_constant_fn`, `test_piecewise_linear_fn`, `test_skyline`, `test_total_lh` |
| [`inference/__tests__/mod.rs`](../../packages/treetime/src/timetree/inference/__tests__/mod.rs)                               | `test_backward_pass`, `test_branch_length_likelihood`, `test_gm_runner`, `test_runner`                                                                                                                                         |
| [`inference/__tests__/test_gm_runner/mod.rs`](../../packages/treetime/src/timetree/inference/__tests__/test_gm_runner/mod.rs) | `test_gm_runner_marginal_dense`, `test_gm_runner_marginal_sparse`, `test_gm_runner_poisson`, `test_gm_runner_pre_optimize`, `test_gm_runner_support`, `test_runner_coalescent`                                                 |
| [`convergence/__tests__/mod.rs`](../../packages/treetime/src/timetree/convergence/__tests__/mod.rs)                           | `test_metrics`, `test_pipeline`                                                                                                                                                                                                |
| [`output/__tests__/mod.rs`](../../packages/treetime/src/commands/timetree/output/__tests__/mod.rs)                                     | `test_auspice`, `test_confidence_combine`, `test_confidence_extract`, `test_confidence_rate`                                                                                                                                   |
| [`optimization/__tests__/mod.rs`](../../packages/treetime/src/timetree/optimization/__tests__/mod.rs)                         | `test_clock_filter`, `test_polytomy`, `test_relaxed_clock`, `test_reroot`                                                                                                                                                      |

---

## Known Test Limitations

1. Golden-master datasets disabled: Multiple datasets disabled due to zero-length branches:
   `dengue_20`, `lassa_L_20`, `mpox_clade_ii_20`, `rsv_a_20`, `tb_20`, `zika_20`

2. Tolerance variations: Root node dominates max diff in Poisson mode tests. Non-root nodes agree within 1e-2.

3. Dense marginal: `ebola_20` disabled due to gap character handling in alphabet.

4. Duplicate tests: [`packages/treetime/src/coalescent/__tests__/test_lineage_dynamics.rs`](../../packages/treetime/src/coalescent/__tests__/test_lineage_dynamics.rs) contains 2 piecewise constant tests that duplicate those in [`packages/treetime-grid/src/__tests__/piecewise_constant_fn.rs`](../../packages/treetime-grid/src/__tests__/piecewise_constant_fn.rs).
