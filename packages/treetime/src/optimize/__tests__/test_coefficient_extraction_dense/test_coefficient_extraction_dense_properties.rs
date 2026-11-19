#[cfg(test)]
mod tests {
  use super::super::test_coefficient_extraction_dense_support::tests::make_dense_seq_dis;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::dense_eval::evaluate_dense_contribution;
  use crate::partition::optimize::dense::get_coefficients;
  use ndarray::{Array1, Axis, array};
  use rstest::rstest;
  use treetime_utils::pretty_assert_abs_diff_eq;

  #[rustfmt::skip]
  #[rstest]
  #[case::zero_partial_overlap(0.0,   array![0.5, 0.3, 0.1, 0.1], array![0.4, 0.4, 0.1, 0.1])]
  #[case::zero_match(          0.0,   array![1.0, 0.0, 0.0, 0.0], array![1.0, 0.0, 0.0, 0.0])]
  #[case::tiny(                0.001, array![0.7, 0.1, 0.1, 0.1], array![0.1, 0.7, 0.1, 0.1])]
  #[case::moderate(            0.1,   array![0.7, 0.1, 0.1, 0.1], array![0.1, 0.7, 0.1, 0.1])]
  #[case::large(               1.0,   array![0.7, 0.1, 0.1, 0.1], array![0.1, 0.7, 0.1, 0.1])]
  #[case::mismatch(            0.1,   array![1.0, 0.0, 0.0, 0.0], array![0.0, 1.0, 0.0, 0.0])]
  #[trace]
  fn test_coefficients_likelihood_matches_jc69_transition_probabilities(
    #[case] branch_length: f64,
    #[case] parent: Array1<f64>,
    #[case] child: Array1<f64>,
  ) {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");
    let overlap = parent.dot(&child);
    let contribution = get_coefficients(
      &make_dense_seq_dis(parent.insert_axis(Axis(0))),
      &make_dense_seq_dis(child.insert_axis(Axis(0))),
      &gtr,
    );

    let metrics = evaluate_dense_contribution(&contribution, branch_length, true).expect("valid branch length");

    let decay = (-4.0 / 3.0 * branch_length).exp();
    let site_lh = 0.25 + decay * (overlap - 0.25);
    let site_lh_derivative = -4.0 / 3.0 * decay * (overlap - 0.25);
    let site_lh_second_derivative = 16.0 / 9.0 * decay * (overlap - 0.25);
    pretty_assert_abs_diff_eq!(site_lh.ln(), metrics.log_lh.value(), epsilon = 1e-12);
    pretty_assert_abs_diff_eq!(site_lh_derivative / site_lh, metrics.derivative, epsilon = 1e-12);
    pretty_assert_abs_diff_eq!(
      site_lh_second_derivative / site_lh - (site_lh_derivative / site_lh).powi(2),
      metrics.second_derivative,
      epsilon = 1e-12
    );
  }
}
