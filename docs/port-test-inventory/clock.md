# Clock Inference Tests

[Back to index](_index.md)

## Summary

| Category                                   | Type                        |
| ------------------------------------------ | --------------------------- |
| [Clock regression](#clock-regression)      | Unit                        |
| [Date constraints](#date-constraints)      | Unit                        |
| [Rerooting](#rerooting)                    | Unit                        |
| [Clock filter](#clock-filter)              | Unit                        |
| [Find best root](#find-best-root)          | Unit                        |
| [Dengue/100 pipeline](#dengue100-pipeline) | Integration + Golden-master |

---

## Clock Regression

**File:** [`test_clock_regression.rs`](../../packages/treetime/src/commands/clock/__tests__/test_clock_regression.rs)

| Test                    | Purpose                                             |
| ----------------------- | --------------------------------------------------- |
| `test_clock_naive_rate` | Naive rate matches WLS regression, weighted variant |

---

## Date Constraints

**File:** [`test_date_constraints.rs`](../../packages/treetime/src/commands/clock/__tests__/test_date_constraints.rs)

| Test                                                       | Purpose                               |
| ---------------------------------------------------------- | ------------------------------------- |
| `test_load_date_constraints_success_three_leaves`          | Three leaves with point dates         |
| `test_load_date_constraints_mixed_leaves`                  | Missing date marks bad branch         |
| `test_load_date_constraints_range`                         | Date range creates Range distribution |
| `test_load_date_constraints_internal_node`                 | Internal node with point date         |
| `test_load_date_constraints_bad_branch_propagation`        | Bad branch on undated leaf            |
| `test_load_date_constraints_all_children_bad`              | Subtree bad branch propagation        |
| `test_load_date_constraints_boundary_exactly_three_leaves` | Minimum three dated leaves            |
| `test_load_date_constraints_date_with_none_value`          | Null date marks bad branch            |
| `test_load_date_constraints_deep_tree_propagation`         | Deep tree bad branch propagation      |
| `test_load_date_constraints_wide_tree`                     | Wide tree with sparse dates           |
| `test_load_date_constraints_mixed_ranges_and_points`       | Mixed point and range dates           |
| `test_load_date_constraints_idempotency`                   | Double load produces same result      |
| `test_load_date_constraints_internal_node_with_range`      | Internal node with range date         |
| `test_load_date_constraints_negative_time`                 | Negative year fractions (BCE dates)   |

---

## Rerooting

**File:** [`test_reroot.rs`](../../packages/treetime/src/commands/clock/__tests__/test_reroot.rs)

| Test                                                                     | Purpose                              |
| ------------------------------------------------------------------------ | ------------------------------------ |
| `test_remove_node_if_trivial_simple`                                     | Remove single-child internal node    |
| `test_reroot_policy_allow_edge_split_false_no_new_nodes`                 | No edge split preserves node count   |
| `test_reroot_policy_remove_old_root_if_trivial_false_preserves_old_root` | Old root preserved when flag off     |
| `test_reroot_policy_default_allows_edge_split`                           | Default policy allows edge splitting |

---

## Find Best Root

**File:** [`test_find_best_root.rs`](../../packages/treetime/src/commands/clock/find_best_root/__tests__/test_find_best_root.rs)

| Test                                                             | Purpose                                                |
| ---------------------------------------------------------------- | ------------------------------------------------------ |
| `test_find_best_root_grid`                                       | Grid search root placement                             |
| `test_find_best_root_grid_with_params`                           | Grid search with custom n_points                       |
| `test_find_best_root_brent`                                      | Brent method root placement                            |
| `test_find_best_root_brent_with_params`                          | Brent method with custom params                        |
| `test_find_best_root_golden_section`                             | Golden section root placement                          |
| `test_find_best_root_golden_section_with_params`                 | Golden section with custom params                      |
| `test_optimization_methods_improve_on_grid`                      | Brent and golden improve on grid                       |
| `test_find_best_root_force_positive_true_rejects_negative_rate`  | Error when all rates negative and force_positive=true  |
| `test_find_best_root_force_positive_false_accepts_negative_rate` | Success when force_positive=false allows negative rate |

---

## Clock Filter

**File:** [`test_clock_filter.rs`](../../packages/treetime/src/commands/clock/__tests__/test_clock_filter.rs)

| Test                                                       | Purpose                                                          |
| ---------------------------------------------------------- | ---------------------------------------------------------------- |
| `test_clock_filter_positive_rate_identifies_outliers`      | IQD filter identifies outliers with positive-rate model          |
| `test_clock_filter_negative_rate_identifies_same_outliers` | Same outliers flagged with negative-rate model (sign-invariance) |

---

## Dengue/100 Pipeline

**File:** [`test_clock_dengue100.rs`](../../packages/treetime/src/commands/clock/__tests__/test_clock_dengue100.rs)

| Test                                                  | Purpose                                                                        |
| ----------------------------------------------------- | ------------------------------------------------------------------------------ |
| `test_dengue100_clock_pipeline_structural_properties` | Assertion-based: rate positive, plausible range, outliers detected, v0 overlap |
| `test_dengue100_clock_pipeline_golden_master`         | Pin v1 output: rate, intercept, R, chisq, outlier set                          |
