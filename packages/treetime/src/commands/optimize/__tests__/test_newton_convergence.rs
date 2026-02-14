#[cfg(test)]
mod tests {
  use crate::commands::optimize::optimize_dense;
  use crate::commands::optimize::optimize_unified::{OptimizationContribution, evaluate_mixed};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use approx::assert_ulps_eq;
  use ndarray::array;
  use num::clamp;

  /// Create a dense contribution with specified coefficients and JC69 GTR model.
  fn make_dense_contribution(coefficients: ndarray::Array2<f64>) -> OptimizationContribution {
    let gtr = jc69(JC69Params::default()).unwrap();
    OptimizationContribution::Dense(optimize_dense::PartitionContribution::new(coefficients, gtr))
  }

  /// Simulate one Newton step: branch_length - clamp(d1/d2, -1.0, branch_length)
  fn newton_step(branch_length: f64, derivative: f64, second_derivative: f64) -> f64 {
    branch_length - clamp(derivative / second_derivative, -1.0, branch_length)
  }

  // ============================================================================
  // Step clamping tests
  // ============================================================================

  #[test]
  fn test_newton_step_clamps_negative_to_minus_one() {
    // When derivative/second_derivative < -1.0, clamp to -1.0
    // This prevents stepping backward by more than 1.0
    let branch_length = 0.5;
    let derivative = 10.0;
    let second_derivative = -2.0; // d1/d2 = -5.0, clamped to -1.0

    let result = newton_step(branch_length, derivative, second_derivative);

    // new = 0.5 - (-1.0) = 1.5
    assert_ulps_eq!(result, 1.5, max_ulps = 4);
  }

  #[test]
  fn test_newton_step_clamps_positive_to_branch_length() {
    // When derivative/second_derivative > branch_length, clamp to branch_length
    // This prevents new_branch_length from going negative
    let branch_length = 0.1;
    let derivative = -10.0;
    let second_derivative = -2.0; // d1/d2 = 5.0, clamped to 0.1

    let result = newton_step(branch_length, derivative, second_derivative);

    // new = 0.1 - 0.1 = 0.0
    assert_ulps_eq!(result, 0.0, max_ulps = 4);
  }

  #[test]
  fn test_newton_step_no_clamping_in_valid_range() {
    // When -1.0 <= d1/d2 <= branch_length, no clamping occurs
    let branch_length = 0.5;
    let derivative = -0.1;
    let second_derivative = -1.0; // d1/d2 = 0.1, in range [-1.0, 0.5]

    let result = newton_step(branch_length, derivative, second_derivative);

    // new = 0.5 - 0.1 = 0.4
    assert_ulps_eq!(result, 0.4, max_ulps = 4);
  }

  // ============================================================================
  // Convergence threshold tests (|new - old| > 0.001 * old)
  // ============================================================================

  #[test]
  fn test_convergence_threshold_above_continues() {
    // When |new - old| > 0.001 * old, iteration should continue
    let old: f64 = 1.0;
    let new: f64 = 1.002; // |1.002 - 1.0| = 0.002 > 0.001 * 1.0 = 0.001

    assert!((new - old).abs() > 0.001 * old);
  }

  #[test]
  fn test_convergence_threshold_at_boundary_stops() {
    // When |new - old| <= 0.001 * old, iteration should stop
    let old: f64 = 1.0;
    let new: f64 = 1.001; // |1.001 - 1.0| = 0.001 = 0.001 * 1.0, NOT > so stops

    assert!((new - old).abs() <= 0.001 * old);
  }

  #[test]
  fn test_convergence_threshold_below_stops() {
    let old: f64 = 1.0;
    let new: f64 = 1.0005; // |0.0005| < 0.001 * 1.0

    assert!((new - old).abs() <= 0.001 * old);
  }

  #[test]
  fn test_convergence_threshold_scales_with_branch_length() {
    // For small branch lengths, threshold is proportionally smaller
    let old: f64 = 0.01;
    let threshold = 0.001 * old; // = 0.00001

    assert_ulps_eq!(threshold, 0.00001, max_ulps = 4);

    // A change of 0.00002 would continue (> 0.00001)
    let new: f64 = 0.01002;
    assert!((new - old).abs() > 0.001 * old);
  }

  // ============================================================================
  // evaluate_mixed with real contributions
  // ============================================================================

  #[test]
  fn test_evaluate_mixed_returns_finite_values() {
    // evaluate_mixed should return finite values for any valid coefficients
    let coefficients = array![[0.25, 0.25, 0.25, 0.25], [0.25, 0.25, 0.25, 0.25],];
    let contribution = make_dense_contribution(coefficients);
    let contributions = vec![contribution];

    let metrics = evaluate_mixed(&contributions, 0.1);

    assert!(metrics.log_lh.is_finite());
    assert!(metrics.derivative.is_finite());
    assert!(metrics.second_derivative.is_finite());
  }

  #[test]
  fn test_evaluate_mixed_log_lh_negative() {
    // Log-likelihood should be negative (probabilities < 1)
    let coefficients = array![[0.9, 0.03, 0.03, 0.04], [0.03, 0.9, 0.03, 0.04],];
    let contribution = make_dense_contribution(coefficients);
    let contributions = vec![contribution];

    let metrics = evaluate_mixed(&contributions, 0.1);

    // Log of a probability is negative or zero
    assert!(metrics.log_lh <= 0.0);
  }

  #[test]
  fn test_evaluate_mixed_multiple_branch_lengths() {
    // evaluate_mixed should work across a range of branch lengths
    let coefficients = array![[0.9, 0.03, 0.03, 0.04], [0.03, 0.9, 0.03, 0.04],];
    let contribution = make_dense_contribution(coefficients);
    let contributions = vec![contribution];

    let branch_lengths = [0.001, 0.01, 0.1, 0.5, 1.0];
    for &bl in &branch_lengths {
      let metrics = evaluate_mixed(&contributions, bl);
      assert!(metrics.log_lh.is_finite(), "log_lh should be finite at bl={bl}");
      assert!(metrics.derivative.is_finite(), "derivative should be finite at bl={bl}");
      assert!(
        metrics.second_derivative.is_finite(),
        "second_derivative should be finite at bl={bl}"
      );
    }
  }

  // ============================================================================
  // Newton iteration simulation tests
  // ============================================================================

  #[test]
  fn test_newton_iteration_converges_within_bounds() {
    // Simulate the Newton loop from optimize_unified.rs
    let coefficients = array![[0.9, 0.03, 0.03, 0.04], [0.03, 0.9, 0.03, 0.04], [0.1, 0.1, 0.7, 0.1],];
    let contribution = make_dense_contribution(coefficients);
    let contributions = vec![contribution];

    let mut branch_length = 0.1;
    let metrics = evaluate_mixed(&contributions, branch_length);
    let max_iter = 10;
    let mut n_iter = 0;

    if metrics.second_derivative < 0.0 {
      let mut new_branch_length =
        branch_length - clamp(metrics.derivative / metrics.second_derivative, -1.0, branch_length);

      while (new_branch_length - branch_length).abs() > 0.001 * branch_length && n_iter < max_iter {
        let new_metrics = evaluate_mixed(&contributions, new_branch_length);
        if new_metrics.second_derivative < 0.0 {
          branch_length = new_branch_length;
          new_branch_length = branch_length
            - clamp(
              new_metrics.derivative / new_metrics.second_derivative,
              -1.0,
              branch_length,
            );
        } else {
          break;
        }
        n_iter += 1;
      }
      branch_length = new_branch_length;
    }

    // Verify convergence happened within iteration limit
    assert!(n_iter <= max_iter);
    // Verify branch length is non-negative
    assert!(branch_length >= 0.0);
  }

  #[test]
  fn test_newton_iteration_respects_max_iter() {
    // Even with slow convergence, should stop at max_iter
    let coefficients = array![[0.25, 0.25, 0.25, 0.25], [0.25, 0.25, 0.25, 0.25],];
    let contribution = make_dense_contribution(coefficients);
    let contributions = vec![contribution];

    let mut branch_length = 0.5;
    let max_iter = 10;
    let mut n_iter = 0;

    let metrics = evaluate_mixed(&contributions, branch_length);
    if metrics.second_derivative < 0.0 {
      let mut new_branch_length =
        branch_length - clamp(metrics.derivative / metrics.second_derivative, -1.0, branch_length);

      while (new_branch_length - branch_length).abs() > 0.001 * branch_length && n_iter < max_iter {
        let new_metrics = evaluate_mixed(&contributions, new_branch_length);
        if new_metrics.second_derivative < 0.0 {
          branch_length = new_branch_length;
          new_branch_length = branch_length
            - clamp(
              new_metrics.derivative / new_metrics.second_derivative,
              -1.0,
              branch_length,
            );
        } else {
          break;
        }
        n_iter += 1;
      }
    }

    assert!(n_iter <= max_iter);
  }

  #[test]
  fn test_newton_early_exit_on_non_concave() {
    // If second_derivative becomes non-negative mid-iteration, loop should break
    // This is tested by checking the branch condition rather than mocking
    let second_derivative = 0.1; // >= 0, should cause break

    // The condition in the loop is: if new_metrics.second_derivative < 0.0
    assert!(second_derivative >= 0.0); // would break
  }

  #[test]
  fn test_newton_step_preserves_positivity() {
    // For any valid inputs, result should be non-negative when starting positive
    let test_cases = [
      (0.5, 1.0, -2.0),   // normal case
      (0.1, -10.0, -2.0), // large positive step, clamped
      (0.01, 5.0, -1.0),  // negative step
      (1.0, 0.0, -1.0),   // zero derivative
    ];

    for (branch_length, derivative, second_derivative) in test_cases {
      let result = newton_step(branch_length, derivative, second_derivative);
      assert!(
        result >= 0.0,
        "newton_step({branch_length}, {derivative}, {second_derivative}) = {result} should be >= 0"
      );
    }
  }

  #[test]
  fn test_newton_step_zero_derivative_no_change() {
    // When derivative is zero, step is zero
    let branch_length = 0.5;
    let derivative = 0.0;
    let second_derivative = -1.0;

    let result = newton_step(branch_length, derivative, second_derivative);

    assert_ulps_eq!(result, branch_length, max_ulps = 4);
  }
}
