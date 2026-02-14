#[cfg(test)]
mod tests {
  use crate::commands::optimize::optimize_dense::{PartitionContribution, evaluate, get_coefficients};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::representation::graph_dense::DenseSeqDis;
  use approx::assert_ulps_eq;
  use ndarray::{Array2, array};

  /// Create a DenseSeqDis from a 2D array of probabilities.
  fn make_dense_seq_dis(dis: Array2<f64>) -> DenseSeqDis {
    DenseSeqDis::new(dis)
  }

  // ==========================================================================
  // Basic coefficient computation tests
  // ==========================================================================

  #[test]
  fn test_get_coefficients_identity_messages() {
    // When both messages are uniform probability, coefficients reflect the
    // eigenvector structure. For JC69 with uniform messages, the coefficient
    // sum equals the transition probability at branch_length=0, which is the
    // product s_a * delta_{ab} * r_b summed over a,b = sum_a (s_a * r_a).
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    // Uniform message: equal probability for all states at each position
    let uniform = array![[0.25, 0.25, 0.25, 0.25]];
    let msg_to_parent = make_dense_seq_dis(uniform.clone());
    let msg_to_child = make_dense_seq_dis(uniform);

    let contribution = get_coefficients(&msg_to_parent, &msg_to_child, &gtr);

    // Coefficients shape should match: 1 position x 4 eigenvalues
    assert_eq!(contribution.coefficients.shape(), &[1, 4]);

    // At branch_length=0, the likelihood is the row sum of coefficients.
    // For uniform messages, this equals sum_a (0.25 * 0.25) = 4 * 0.0625 = 0.25
    let row_sum = contribution.coefficients.sum_axis(ndarray::Axis(1));
    assert_ulps_eq!(row_sum[0], 0.25, max_ulps = 10);

    // Verify via evaluate function
    let metrics = evaluate(&[contribution], 0.0);
    assert_ulps_eq!(metrics.log_lh, 0.25_f64.ln(), max_ulps = 100);
  }

  #[test]
  fn test_get_coefficients_certain_state_parent() {
    // Parent message is certain (state A), child is uniform
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let parent = array![[1.0, 0.0, 0.0, 0.0]]; // certain A
    let child = array![[0.25, 0.25, 0.25, 0.25]]; // uniform
    let msg_to_parent = make_dense_seq_dis(parent);
    let msg_to_child = make_dense_seq_dis(child);

    let contribution = get_coefficients(&msg_to_parent, &msg_to_child, &gtr);

    // Row sum should still sum to product of individual sums
    let row_sum = contribution.coefficients.sum_axis(ndarray::Axis(1))[0];
    assert_ulps_eq!(row_sum, 0.25, max_ulps = 10);
  }

  #[test]
  fn test_get_coefficients_certain_state_child() {
    // Parent is uniform, child message is certain (state A)
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let parent = array![[0.25, 0.25, 0.25, 0.25]]; // uniform
    let child = array![[1.0, 0.0, 0.0, 0.0]]; // certain A
    let msg_to_parent = make_dense_seq_dis(parent);
    let msg_to_child = make_dense_seq_dis(child);

    let contribution = get_coefficients(&msg_to_parent, &msg_to_child, &gtr);

    // Row sum at branch_length=0 equals sum_a (child_a * parent_a) = 1.0 * 0.25 = 0.25
    let row_sum = contribution.coefficients.sum_axis(ndarray::Axis(1))[0];
    assert_ulps_eq!(row_sum, 0.25, max_ulps = 10);
  }

  #[test]
  fn test_get_coefficients_matching_certain_states() {
    // Both parent and child are certain of same state A
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let certain_a = array![[1.0, 0.0, 0.0, 0.0]];
    let msg_to_parent = make_dense_seq_dis(certain_a.clone());
    let msg_to_child = make_dense_seq_dis(certain_a);

    let contribution = get_coefficients(&msg_to_parent, &msg_to_child, &gtr);

    // At branch_length=0, likelihood should be high (states match)
    let metrics = evaluate(&[contribution], 0.0);
    assert!(metrics.log_lh > -1.0, "log-LH should be high for matching states");
  }

  #[test]
  fn test_get_coefficients_mismatched_certain_states() {
    // Parent certain A, child certain C - require substitution
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let parent = array![[1.0, 0.0, 0.0, 0.0]]; // certain A
    let child = array![[0.0, 1.0, 0.0, 0.0]]; // certain C
    let msg_to_parent = make_dense_seq_dis(parent);
    let msg_to_child = make_dense_seq_dis(child);

    let contribution = get_coefficients(&msg_to_parent, &msg_to_child, &gtr);

    // At branch_length=0, likelihood should be zero (mismatch with no time for substitution)
    let metrics = evaluate(&[contribution], 0.0);
    assert!(
      metrics.log_lh < -10.0 || metrics.log_lh == f64::NEG_INFINITY,
      "log-LH should be very low for mismatched states at zero branch length"
    );
  }

  // ==========================================================================
  // Multiple positions tests
  // ==========================================================================

  #[test]
  fn test_get_coefficients_multiple_positions() {
    // Test with multiple alignment positions
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let parent = array![
      [1.0, 0.0, 0.0, 0.0],     // position 0: certain A
      [0.0, 1.0, 0.0, 0.0],     // position 1: certain C
      [0.25, 0.25, 0.25, 0.25]  // position 2: uniform
    ];
    let child = array![
      [1.0, 0.0, 0.0, 0.0],     // position 0: certain A (match)
      [0.0, 1.0, 0.0, 0.0],     // position 1: certain C (match)
      [0.25, 0.25, 0.25, 0.25]  // position 2: uniform
    ];
    let msg_to_parent = make_dense_seq_dis(parent);
    let msg_to_child = make_dense_seq_dis(child);

    let contribution = get_coefficients(&msg_to_parent, &msg_to_child, &gtr);

    // Should have 3 rows of coefficients
    assert_eq!(contribution.coefficients.nrows(), 3);
    assert_eq!(contribution.coefficients.ncols(), 4); // 4 eigenvalues for JC69
  }

  #[test]
  fn test_get_coefficients_row_independence() {
    // Each position's coefficients should be computed independently
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    // Single position with certain A
    let single_parent = array![[1.0, 0.0, 0.0, 0.0]];
    let single_child = array![[0.25, 0.25, 0.25, 0.25]];
    let single_contribution = get_coefficients(
      &make_dense_seq_dis(single_parent),
      &make_dense_seq_dis(single_child),
      &gtr,
    );

    // Multiple positions, first matches the single case
    let multi_parent = array![
      [1.0, 0.0, 0.0, 0.0],
      [0.0, 0.0, 1.0, 0.0] // different state
    ];
    let multi_child = array![
      [0.25, 0.25, 0.25, 0.25],
      [0.5, 0.5, 0.0, 0.0] // different distribution
    ];
    let multi_contribution = get_coefficients(
      &make_dense_seq_dis(multi_parent),
      &make_dense_seq_dis(multi_child),
      &gtr,
    );

    // First row should match single contribution
    for i in 0..4 {
      assert_ulps_eq!(
        single_contribution.coefficients[[0, i]],
        multi_contribution.coefficients[[0, i]],
        max_ulps = 10
      );
    }
  }

  // ==========================================================================
  // Mathematical property tests
  // ==========================================================================

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

  // ==========================================================================
  // GTR matrix property tests
  // ==========================================================================

  #[test]
  fn test_coefficients_use_eigenvector_decomposition() {
    // Verify that coefficients are computed via eigenvector decomposition:
    // k_c = (child · v)_c * (parent · v_inv^T)_c
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let parent = array![[0.4, 0.3, 0.2, 0.1]];
    let child = array![[0.1, 0.2, 0.3, 0.4]];
    let msg_to_parent = make_dense_seq_dis(parent.clone());
    let msg_to_child = make_dense_seq_dis(child.clone());

    let contribution = get_coefficients(&msg_to_parent, &msg_to_child, &gtr);

    // Manually compute expected coefficients
    let child_v = child.dot(&gtr.v);
    let parent_v_inv_t = parent.dot(&gtr.v_inv.t());
    let expected = child_v * parent_v_inv_t;

    for i in 0..4 {
      assert_ulps_eq!(contribution.coefficients[[0, i]], expected[[0, i]], max_ulps = 10);
    }
  }

  #[test]
  fn test_coefficients_preserve_gtr_reference() {
    // PartitionContribution should store the GTR for later evaluation
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let parent = array![[0.25, 0.25, 0.25, 0.25]];
    let child = array![[0.25, 0.25, 0.25, 0.25]];
    let msg_to_parent = make_dense_seq_dis(parent);
    let msg_to_child = make_dense_seq_dis(child);

    let contribution = get_coefficients(&msg_to_parent, &msg_to_child, &gtr);

    // Stored GTR should have same eigenvalues
    for i in 0..4 {
      assert_ulps_eq!(contribution.gtr.eigvals[i], gtr.eigvals[i], max_ulps = 10);
    }
  }

  // ==========================================================================
  // PartitionContribution::new tests
  // ==========================================================================

  #[test]
  fn test_partition_contribution_new() {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");
    let coefficients = array![[0.1, 0.2, 0.3, 0.4], [0.4, 0.3, 0.2, 0.1]];

    let contribution = PartitionContribution::new(coefficients.clone(), gtr);

    assert_eq!(contribution.coefficients.shape(), coefficients.shape());
    for i in 0..2 {
      for j in 0..4 {
        assert_ulps_eq!(contribution.coefficients[[i, j]], coefficients[[i, j]], max_ulps = 10);
      }
    }
  }
}
