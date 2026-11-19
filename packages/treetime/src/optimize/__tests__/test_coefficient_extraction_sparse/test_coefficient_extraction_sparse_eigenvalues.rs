#[cfg(test)]
mod tests {
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::sparse_eval::evaluate_sparse_contribution;
  use crate::partition::optimize::sparse::{PartitionContribution, SiteContribution};
  use ndarray::array;
  use rstest::rstest;
  use treetime_utils::pretty_assert_abs_diff_eq;

  #[rustfmt::skip]
  #[rstest]
  #[case::short(0.01)]
  #[case::long( 1.0)]
  #[trace]
  fn test_sparse_evaluation_matches_jc69_identity_probability(#[case] branch_length: f64) {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");
    let site = SiteContribution {
      multiplicity: 2.0,
      coefficients: array![0.25, 0.25, 0.25, 0.25],
    };
    let contribution = PartitionContribution::new(vec![site], &gtr);

    let metrics = evaluate_sparse_contribution(&contribution, branch_length, false).expect("valid branch length");

    let jc69_identity_probability = 0.25 + 0.75 * (-4.0 / 3.0 * branch_length).exp();
    pretty_assert_abs_diff_eq!(
      2.0 * jc69_identity_probability.ln(),
      metrics.log_lh.value(),
      epsilon = 1e-12
    );
  }
}
