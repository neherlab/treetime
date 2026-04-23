#[cfg(test)]
mod tests {
  use crate::commands::clock::clock_set::ClockSet;
  use ndarray::array;
  use treetime_utils::pretty_assert_abs_diff_eq;

  #[test]
  fn test_clock_set_cov_is_hessian_inverse() {
    let cs = ClockSet::leaf_contribution_to_parent(Some(2020.5), 0.03, 0.001)
      + ClockSet::leaf_contribution_to_parent(Some(2021.0), 0.05, 0.002)
      + ClockSet::leaf_contribution_to_parent(Some(2019.0), 0.01, 0.001);

    let h = cs.hessian();
    let c = cs.cov();

    // H * H^{-1} = I
    let product = h.dot(&c);
    let identity = array![[1.0, 0.0], [0.0, 1.0]];
    pretty_assert_abs_diff_eq!(product, identity, epsilon = 1e-6);
  }

  #[test]
  fn test_clock_set_cov_off_diagonal_symmetry() {
    let cs = ClockSet::leaf_contribution_to_parent(Some(2020.5), 0.03, 0.001)
      + ClockSet::leaf_contribution_to_parent(Some(2021.0), 0.05, 0.002);

    let c = cs.cov();
    pretty_assert_abs_diff_eq!(c[[0, 1]], c[[1, 0]], epsilon = 1e-15);
  }

  #[test]
  fn test_clock_set_cov_diagonal_positive() {
    let cs = ClockSet::leaf_contribution_to_parent(Some(2020.5), 0.03, 0.001)
      + ClockSet::leaf_contribution_to_parent(Some(2021.0), 0.05, 0.002);

    let c = cs.cov();
    assert!(c[[0, 0]] > 0.0, "rate variance must be positive");
    assert!(c[[1, 1]] > 0.0, "intercept variance must be positive");
  }

  #[test]
  fn test_clock_set_cov_off_diagonal_uses_t_sum_not_dt_sum() {
    // Construct a ClockSet where t_sum != dt_sum to distinguish the two
    let cs = ClockSet::leaf_contribution_to_parent(Some(2020.0), 0.05, 0.001)
      + ClockSet::leaf_contribution_to_parent(Some(2022.0), 0.01, 0.002);

    assert!(
      (cs.t_sum() - cs.dt_sum()).abs() > 1.0,
      "test requires t_sum != dt_sum to be discriminating"
    );

    // 2x2 inverse: [[a,b],[c,d]]^{-1} = 1/det * [[d,-b],[-c,a]]
    // H = [[tsq_sum, t_sum], [t_sum, norm]]
    // H^{-1}[0,1] = -t_sum / det, NOT -dt_sum / det
    let det = cs.determinant();
    let expected_off_diag = -cs.t_sum() / det;
    let c = cs.cov();
    pretty_assert_abs_diff_eq!(c[[0, 1]], expected_off_diag, epsilon = 1e-10);
  }
}
