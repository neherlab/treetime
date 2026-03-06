use crate::commands::optimize::optimize_dense::{evaluate, get_coefficients};
use crate::gtr::get_gtr::{JC69Params, jc69};
use approx::assert_ulps_eq;
use ndarray::{Axis, array};

use super::test_coefficient_extraction_dense_support::make_dense_seq_dis;

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
  let row_sum = contribution.coefficients.sum_axis(Axis(1));
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
  let row_sum = contribution.coefficients.sum_axis(Axis(1))[0];
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
  let row_sum = contribution.coefficients.sum_axis(Axis(1))[0];
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
