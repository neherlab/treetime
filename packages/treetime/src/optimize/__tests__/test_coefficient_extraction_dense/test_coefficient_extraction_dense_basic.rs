#[cfg(test)]
mod tests {
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::dense_eval::evaluate_dense_contribution;
  use crate::partition::optimize::dense::get_coefficients;
  use crate::pretty_assert_ulps_eq;
  use ndarray::{Axis, array};

  use super::super::test_coefficient_extraction_dense_support::tests::make_dense_seq_dis;

  #[test]
  fn test_get_coefficients_identity_messages() {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let uniform = array![[0.25, 0.25, 0.25, 0.25]];
    let msg_to_parent = make_dense_seq_dis(uniform.clone());
    let msg_to_child = make_dense_seq_dis(uniform);

    let contribution = get_coefficients(&msg_to_parent, &msg_to_child, &gtr);

    assert_eq!(&[1, 4], contribution.coefficients.shape());

    let row_sum = contribution.coefficients.sum_axis(Axis(1));
    pretty_assert_ulps_eq!(row_sum[0], 0.25, max_ulps = 10);

    let metrics = evaluate_dense_contribution(&contribution, 0.0, true).expect("valid branch length");
    pretty_assert_ulps_eq!(metrics.log_lh.value(), 0.25_f64.ln(), max_ulps = 100);
  }

  #[test]
  fn test_get_coefficients_certain_state_parent() {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let parent = array![[1.0, 0.0, 0.0, 0.0]];
    let child = array![[0.25, 0.25, 0.25, 0.25]];
    let msg_to_parent = make_dense_seq_dis(parent);
    let msg_to_child = make_dense_seq_dis(child);

    let contribution = get_coefficients(&msg_to_parent, &msg_to_child, &gtr);

    let row_sum = contribution.coefficients.sum_axis(Axis(1))[0];
    pretty_assert_ulps_eq!(row_sum, 0.25, max_ulps = 10);
  }

  #[test]
  fn test_get_coefficients_certain_state_child() {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let parent = array![[0.25, 0.25, 0.25, 0.25]];
    let child = array![[1.0, 0.0, 0.0, 0.0]];
    let msg_to_parent = make_dense_seq_dis(parent);
    let msg_to_child = make_dense_seq_dis(child);

    let contribution = get_coefficients(&msg_to_parent, &msg_to_child, &gtr);

    let row_sum = contribution.coefficients.sum_axis(Axis(1))[0];
    pretty_assert_ulps_eq!(row_sum, 0.25, max_ulps = 10);
  }

  #[test]
  fn test_get_coefficients_matching_certain_states() {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let certain_a = array![[1.0, 0.0, 0.0, 0.0]];
    let msg_to_parent = make_dense_seq_dis(certain_a.clone());
    let msg_to_child = make_dense_seq_dis(certain_a);

    let contribution = get_coefficients(&msg_to_parent, &msg_to_child, &gtr);

    let metrics = evaluate_dense_contribution(&contribution, 0.0, true).expect("valid branch length");
    assert!(
      metrics.log_lh.value() > -1.0,
      "log-LH should be high for matching states"
    );
  }

  #[test]
  fn test_get_coefficients_mismatched_certain_states() {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let parent = array![[1.0, 0.0, 0.0, 0.0]];
    let child = array![[0.0, 1.0, 0.0, 0.0]];
    let msg_to_parent = make_dense_seq_dis(parent);
    let msg_to_child = make_dense_seq_dis(child);

    let contribution = get_coefficients(&msg_to_parent, &msg_to_child, &gtr);

    let metrics = evaluate_dense_contribution(&contribution, 0.0, true).expect("valid branch length");
    assert!(
      metrics.log_lh.value() < -10.0 || metrics.log_lh.value() == f64::NEG_INFINITY,
      "log-LH should be very low for mismatched states at zero branch length"
    );
  }
}
