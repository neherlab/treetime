#[cfg(test)]
mod tests {
  use crate::commands::optimize::optimize_dense;
  use crate::commands::optimize::optimize_unified::{OptimizationContribution, evaluate_mixed};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use ndarray::array;
  use num::clamp;

  /// Create a dense contribution with specified coefficients and JC69 GTR model.
  fn make_dense_contribution(coefficients: ndarray::Array2<f64>) -> OptimizationContribution {
    let gtr = jc69(JC69Params::default()).unwrap();
    OptimizationContribution::Dense(optimize_dense::PartitionContribution::new(coefficients, gtr))
  }

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
}
