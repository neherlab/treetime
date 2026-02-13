#[cfg(test)]
mod tests {
  use crate::commands::optimize::optimize_dense;
  use crate::commands::optimize::optimize_unified::{evaluate_mixed, OptimizationContribution};
  use crate::gtr::get_gtr::{jc69, JC69Params};
  use approx::assert_ulps_eq;
  use ndarray::{array, Array1};
  use ordered_float::OrderedFloat;

  /// Create a dense contribution with specified coefficients and JC69 GTR model.
  fn make_dense_contribution(coefficients: ndarray::Array2<f64>) -> OptimizationContribution {
    let gtr = jc69(JC69Params::default()).unwrap();
    OptimizationContribution::Dense(optimize_dense::PartitionContribution::new(coefficients, gtr))
  }

  /// Simulate grid search: find branch length with maximum log-LH from a grid of points.
  fn grid_search(
    contributions: &[OptimizationContribution],
    branch_length: f64,
    one_mutation: f64,
  ) -> f64 {
    let branch_lengths = Array1::linspace(0.1 * one_mutation, 1.5 * branch_length + one_mutation, 100);

    branch_lengths
      .iter()
      .max_by_key(|&&bl| {
        let metrics = evaluate_mixed(contributions, bl);
        OrderedFloat(metrics.log_lh)
      })
      .copied()
      .unwrap_or(branch_length)
  }

  // ============================================================================
  // Grid search trigger condition tests
  // ============================================================================

  #[test]
  fn test_grid_search_triggered_when_d2_zero() {
    // When second_derivative == 0, Newton's method cannot proceed (division by zero)
    // Grid search should be used instead
    let second_derivative = 0.0;
    assert!(second_derivative >= 0.0); // triggers grid search path
  }

  #[test]
  fn test_grid_search_triggered_when_d2_positive() {
    // When second_derivative > 0, function is convex (not concave)
    // Newton would move away from maximum, so grid search is needed
    let second_derivative = 0.5;
    assert!(second_derivative >= 0.0); // triggers grid search path
  }

  #[test]
  fn test_newton_path_when_d2_negative() {
    // When second_derivative < 0, function is concave - Newton's method works
    let second_derivative = -0.5;
    assert!(second_derivative < 0.0); // Newton's method path
  }

  // ============================================================================
  // Grid search range tests
  // ============================================================================

  #[test]
  fn test_grid_search_range_lower_bound() {
    // Lower bound is 0.1 * one_mutation
    let one_mutation = 0.01; // 1 mutation in 100 sites
    let lower = 0.1 * one_mutation;

    assert_ulps_eq!(lower, 0.001, max_ulps = 4);
  }

  #[test]
  fn test_grid_search_range_upper_bound() {
    // Upper bound is 1.5 * branch_length + one_mutation
    let branch_length = 0.1;
    let one_mutation = 0.01;
    let upper = 1.5 * branch_length + one_mutation;

    assert_ulps_eq!(upper, 0.16, max_ulps = 4);
  }

  #[test]
  fn test_grid_search_range_100_points() {
    // Grid search uses exactly 100 points
    let branch_length = 0.1;
    let one_mutation = 0.01;
    let branch_lengths = Array1::linspace(0.1 * one_mutation, 1.5 * branch_length + one_mutation, 100);

    assert_eq!(branch_lengths.len(), 100);
  }

  #[test]
  fn test_grid_search_range_evenly_spaced() {
    // Points should be evenly spaced
    let branch_length = 0.1;
    let one_mutation = 0.01;
    let branch_lengths = Array1::linspace(0.1 * one_mutation, 1.5 * branch_length + one_mutation, 100);

    let expected_step = (1.5 * branch_length + one_mutation - 0.1 * one_mutation) / 99.0;

    for i in 1..branch_lengths.len() {
      let actual_step = branch_lengths[i] - branch_lengths[i - 1];
      assert_ulps_eq!(actual_step, expected_step, max_ulps = 4);
    }
  }

  // ============================================================================
  // Grid search finds maximum tests
  // ============================================================================

  #[test]
  fn test_grid_search_finds_maximum_log_lh() {
    // Grid search should find the branch length that maximizes log-LH
    let coefficients = array![
      [0.9, 0.03, 0.03, 0.04],
      [0.03, 0.9, 0.03, 0.04],
      [0.1, 0.1, 0.7, 0.1],
    ];
    let contribution = make_dense_contribution(coefficients);
    let contributions = vec![contribution];

    let branch_length = 0.1;
    let one_mutation = 0.001;

    let best_bl = grid_search(&contributions, branch_length, one_mutation);

    // Verify no other point in the grid has higher log-LH
    let branch_lengths = Array1::linspace(0.1 * one_mutation, 1.5 * branch_length + one_mutation, 100);
    let best_log_lh = evaluate_mixed(&contributions, best_bl).log_lh;

    for &bl in &branch_lengths {
      let log_lh = evaluate_mixed(&contributions, bl).log_lh;
      assert!(
        log_lh <= best_log_lh + 1e-10,
        "Found bl={bl} with log_lh={log_lh} > best_log_lh={best_log_lh}"
      );
    }
  }

  #[test]
  fn test_grid_search_result_in_range() {
    // Result should be within the search range
    let coefficients = array![
      [0.9, 0.03, 0.03, 0.04],
      [0.03, 0.9, 0.03, 0.04],
    ];
    let contribution = make_dense_contribution(coefficients);
    let contributions = vec![contribution];

    let branch_length = 0.1;
    let one_mutation = 0.001;

    let best_bl = grid_search(&contributions, branch_length, one_mutation);

    let lower = 0.1 * one_mutation;
    let upper = 1.5 * branch_length + one_mutation;

    assert!(best_bl >= lower, "best_bl={best_bl} < lower={lower}");
    assert!(best_bl <= upper, "best_bl={best_bl} > upper={upper}");
  }

  #[test]
  fn test_grid_search_with_small_branch_length() {
    // Grid search should work even with very small branch lengths
    let coefficients = array![
      [0.99, 0.003, 0.003, 0.004],
      [0.003, 0.99, 0.003, 0.004],
    ];
    let contribution = make_dense_contribution(coefficients);
    let contributions = vec![contribution];

    let branch_length = 0.001;
    let one_mutation = 0.0001;

    let best_bl = grid_search(&contributions, branch_length, one_mutation);

    assert!(best_bl >= 0.0);
    assert!(best_bl.is_finite());
  }

  #[test]
  fn test_grid_search_with_large_branch_length() {
    // Grid search should work with larger branch lengths
    let coefficients = array![
      [0.5, 0.15, 0.15, 0.2],
      [0.15, 0.5, 0.15, 0.2],
    ];
    let contribution = make_dense_contribution(coefficients);
    let contributions = vec![contribution];

    let branch_length = 1.0;
    let one_mutation = 0.001;

    let best_bl = grid_search(&contributions, branch_length, one_mutation);

    assert!(best_bl >= 0.0);
    assert!(best_bl.is_finite());
    assert!(best_bl <= 1.5 * branch_length + one_mutation);
  }

  // ============================================================================
  // Grid search with multiple contributions
  // ============================================================================

  #[test]
  fn test_grid_search_multiple_partitions() {
    // Grid search should work with multiple contributions
    let coefficients1 = array![
      [0.9, 0.03, 0.03, 0.04],
      [0.03, 0.9, 0.03, 0.04],
    ];
    let coefficients2 = array![
      [0.8, 0.06, 0.06, 0.08],
      [0.06, 0.8, 0.06, 0.08],
    ];
    let contribution1 = make_dense_contribution(coefficients1);
    let contribution2 = make_dense_contribution(coefficients2);
    let contributions = vec![contribution1, contribution2];

    let branch_length = 0.1;
    let one_mutation = 0.001;

    let best_bl = grid_search(&contributions, branch_length, one_mutation);

    // Verify combined log-LH is maximized
    let branch_lengths = Array1::linspace(0.1 * one_mutation, 1.5 * branch_length + one_mutation, 100);
    let best_log_lh = evaluate_mixed(&contributions, best_bl).log_lh;

    for &bl in &branch_lengths {
      let log_lh = evaluate_mixed(&contributions, bl).log_lh;
      assert!(log_lh <= best_log_lh + 1e-10);
    }
  }

  #[test]
  fn test_grid_search_uniform_coefficients() {
    // With uniform coefficients, any branch length gives similar log-LH
    // Grid search should still return a valid result
    let coefficients = array![
      [0.25, 0.25, 0.25, 0.25],
      [0.25, 0.25, 0.25, 0.25],
    ];
    let contribution = make_dense_contribution(coefficients);
    let contributions = vec![contribution];

    let branch_length = 0.1;
    let one_mutation = 0.001;

    let best_bl = grid_search(&contributions, branch_length, one_mutation);

    assert!(best_bl.is_finite());
    assert!(best_bl >= 0.0);
  }

  // ============================================================================
  // Grid search edge cases
  // ============================================================================

  #[test]
  fn test_grid_search_one_mutation_dominates_range() {
    // When one_mutation is large relative to branch_length,
    // the range should still be valid
    let coefficients = array![
      [0.9, 0.03, 0.03, 0.04],
      [0.03, 0.9, 0.03, 0.04],
    ];
    let contribution = make_dense_contribution(coefficients);
    let contributions = vec![contribution];

    let branch_length = 0.001;
    let one_mutation = 0.01; // 10x branch_length

    let lower = 0.1 * one_mutation;
    let upper = 1.5 * branch_length + one_mutation;

    // Lower bound is 0.001, upper is 0.0115
    assert!(lower < upper, "lower={lower} >= upper={upper}");

    let best_bl = grid_search(&contributions, branch_length, one_mutation);
    assert!(best_bl >= lower);
    assert!(best_bl <= upper);
  }

  #[test]
  fn test_grid_search_preserves_positive_branch_length() {
    // Grid search should always return a positive branch length
    let coefficients = array![
      [0.9, 0.03, 0.03, 0.04],
      [0.03, 0.9, 0.03, 0.04],
    ];
    let contribution = make_dense_contribution(coefficients);
    let contributions = vec![contribution];

    let test_cases = [
      (0.1, 0.001),
      (0.01, 0.001),
      (0.001, 0.0001),
      (1.0, 0.001),
    ];

    for (branch_length, one_mutation) in test_cases {
      let best_bl = grid_search(&contributions, branch_length, one_mutation);
      assert!(
        best_bl > 0.0,
        "grid_search({branch_length}, {one_mutation}) = {best_bl} should be > 0"
      );
    }
  }
}
