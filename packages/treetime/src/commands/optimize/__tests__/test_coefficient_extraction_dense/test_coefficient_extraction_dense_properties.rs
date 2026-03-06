#[cfg(test)]
mod tests {
  use crate::commands::optimize::optimize_dense::{evaluate, get_coefficients};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use ndarray::array;

  use super::super::test_coefficient_extraction_dense_support::tests::make_dense_seq_dis;

  #[test]
  fn test_coefficients_produce_valid_likelihood_at_zero() {
    // At branch_length=0, exp(λt) = 1 for all eigenvalues
    // So likelihood = sum of coefficients per row
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let parent = array![[0.5, 0.3, 0.1, 0.1]];
    let child = array![[0.4, 0.4, 0.1, 0.1]];
    let msg_to_parent = make_dense_seq_dis(parent);
    let msg_to_child = make_dense_seq_dis(child);

    let contribution = get_coefficients(&msg_to_parent, &msg_to_child, &gtr);

    // Evaluate at branch_length=0
    let metrics = evaluate(&[contribution], 0.0);

    // Log-LH should be finite for valid probability distributions
    assert!(metrics.log_lh.is_finite(), "log-LH should be finite");
  }

  #[test]
  fn test_coefficients_likelihood_decreases_for_mismatch_at_zero() {
    // When states mismatch, likelihood at zero branch length should be lower
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    // Matching states
    let match_parent = array![[1.0, 0.0, 0.0, 0.0]];
    let match_child = array![[1.0, 0.0, 0.0, 0.0]];
    let match_contribution = get_coefficients(
      &make_dense_seq_dis(match_parent),
      &make_dense_seq_dis(match_child),
      &gtr,
    );
    let match_metrics = evaluate(&[match_contribution], 0.0);

    // Mismatching states
    let mismatch_parent = array![[1.0, 0.0, 0.0, 0.0]];
    let mismatch_child = array![[0.0, 1.0, 0.0, 0.0]];
    let mismatch_contribution = get_coefficients(
      &make_dense_seq_dis(mismatch_parent),
      &make_dense_seq_dis(mismatch_child),
      &gtr,
    );
    let mismatch_metrics = evaluate(&[mismatch_contribution], 0.0);

    assert!(
      match_metrics.log_lh > mismatch_metrics.log_lh,
      "matching states should have higher log-LH at zero branch length"
    );
  }

  #[test]
  fn test_coefficients_support_positive_branch_length() {
    // Coefficients should produce valid metrics at positive branch lengths
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let parent = array![[0.7, 0.1, 0.1, 0.1]];
    let child = array![[0.1, 0.7, 0.1, 0.1]];
    let msg_to_parent = make_dense_seq_dis(parent);
    let msg_to_child = make_dense_seq_dis(child);

    let contribution = get_coefficients(&msg_to_parent, &msg_to_child, &gtr);
    let contributions = [contribution];

    // Test at various branch lengths
    for &branch_length in &[0.001, 0.01, 0.1, 1.0] {
      let metrics = evaluate(&contributions, branch_length);
      assert!(
        metrics.log_lh.is_finite(),
        "log-LH should be finite at branch_length={branch_length}"
      );
      assert!(
        metrics.derivative.is_finite(),
        "derivative should be finite at branch_length={branch_length}"
      );
      assert!(
        metrics.second_derivative.is_finite(),
        "second_derivative should be finite at branch_length={branch_length}"
      );
    }
  }
}
