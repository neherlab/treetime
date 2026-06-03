#[cfg(test)]
mod tests {
  use crate::partition::marginal_core::normalize_1d_inplace;
  use ndarray::{Array1, array};
  use rstest::rstest;
  use treetime_utils::pretty_assert_abs_diff_eq;

  // Valid (positive, finite) norm: distribution is normalized to sum 1 and the
  // contribution is weight * ln(norm).
  #[rustfmt::skip]
  #[rstest]
  #[case::unweighted(    array![1.0, 3.0],      1.0, array![0.25, 0.75],      4.0_f64.ln())]
  #[case::weighted(      array![1.0, 1.0],      5.0, array![0.5,  0.5],  5.0 * 2.0_f64.ln())]
  #[case::already_norm(  array![0.25, 0.75],    1.0, array![0.25, 0.75],      1.0_f64.ln())]
  #[case::three_states(  array![2.0, 2.0, 4.0], 1.0, array![0.25, 0.25, 0.5], 8.0_f64.ln())]
  #[trace]
  fn test_marginal_core_normalize_1d_inplace_valid(
    #[case] dis: Array1<f64>,
    #[case] weight: f64,
    #[case] expected_dis: Array1<f64>,
    #[case] expected_ll: f64,
  ) {
    let mut dis = dis;
    let ll = normalize_1d_inplace(&mut dis, weight);
    pretty_assert_abs_diff_eq!(expected_dis, dis, epsilon = 1e-12);
    pretty_assert_abs_diff_eq!(expected_ll, ll, epsilon = 1e-12);
  }

  // Non-positive or non-finite norm: fall back to a uniform distribution and
  // contribute NEG_INFINITY (unweighted, regardless of the weight argument).
  #[rustfmt::skip]
  #[rstest]
  #[case::zero_norm(     array![0.0, 0.0, 0.0],          3.0)]
  #[case::infinite_norm( array![f64::INFINITY, 1.0],     2.0)]
  #[case::nan_norm(      array![f64::NAN, 1.0],          1.0)]
  #[trace]
  fn test_marginal_core_normalize_1d_inplace_fallback(#[case] dis: Array1<f64>, #[case] weight: f64) {
    let mut dis = dis;
    let n = dis.len();
    let ll = normalize_1d_inplace(&mut dis, weight);
    let expected_uniform = Array1::from_elem(n, 1.0 / n as f64);
    pretty_assert_abs_diff_eq!(expected_uniform, dis, epsilon = 1e-12);
    assert_eq!(f64::NEG_INFINITY, ll);
  }
}
