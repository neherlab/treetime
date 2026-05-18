#[cfg(test)]
mod tests {
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::pretty_assert_ulps_eq;
  use crate::partition::optimize_dense::get_coefficients;
  use ndarray::array;

  use super::super::test_coefficient_extraction_dense_support::tests::make_dense_seq_dis;

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
      pretty_assert_ulps_eq!(
        single_contribution.coefficients[[0, i]],
        multi_contribution.coefficients[[0, i]],
        max_ulps = 10
      );
    }
  }
}
