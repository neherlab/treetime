#[cfg(test)]
mod tests {
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::optimize::dense::get_coefficients;
  use ndarray::array;
  use treetime_utils::pretty_assert_ulps_eq;

  use super::super::test_coefficient_extraction_dense_support::tests::make_dense_seq_dis;

  #[test]
  fn test_get_coefficients_multiple_positions() {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let parent = array![[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.25, 0.25, 0.25, 0.25]];
    let child = array![[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.25, 0.25, 0.25, 0.25]];
    let msg_to_parent = make_dense_seq_dis(parent);
    let msg_to_child = make_dense_seq_dis(child);

    let contribution = get_coefficients(&msg_to_parent, &msg_to_child, &gtr);

    assert_eq!(3, contribution.coefficients.nrows());
    assert_eq!(4, contribution.coefficients.ncols());
  }

  #[test]
  fn test_get_coefficients_row_independence() {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let single_parent = array![[1.0, 0.0, 0.0, 0.0]];
    let single_child = array![[0.25, 0.25, 0.25, 0.25]];
    let single_contribution = get_coefficients(
      &make_dense_seq_dis(single_parent),
      &make_dense_seq_dis(single_child),
      &gtr,
    );

    let multi_parent = array![[1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]];
    let multi_child = array![[0.25, 0.25, 0.25, 0.25], [0.5, 0.5, 0.0, 0.0]];
    let multi_contribution = get_coefficients(
      &make_dense_seq_dis(multi_parent),
      &make_dense_seq_dis(multi_child),
      &gtr,
    );

    pretty_assert_ulps_eq!(
      single_contribution.coefficients.row(0).to_owned(),
      multi_contribution.coefficients.row(0).to_owned(),
      max_ulps = 10
    );
  }
}
