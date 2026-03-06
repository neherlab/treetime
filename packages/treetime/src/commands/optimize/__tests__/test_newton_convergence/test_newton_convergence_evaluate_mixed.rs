use crate::commands::optimize::optimize_unified::evaluate_mixed;
use ndarray::array;

use super::test_newton_convergence_support::make_dense_contribution;

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
