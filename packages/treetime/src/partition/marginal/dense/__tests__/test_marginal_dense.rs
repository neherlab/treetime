#[cfg(test)]
mod tests {
  use crate::partition::marginal::shared::normalize::{normalize_inplace, normalize_sites_from_log};
  use crate::pretty_assert_neg_inf;
  use helpers::assert_valid_rows;
  use ndarray::{Array2, array};
  use rstest::rstest;
  use treetime_utils::pretty_assert_abs_diff_eq;

  #[test]
  fn test_normalize_sites_from_log_gives_uniform_row_and_neg_inf_norm_for_impossible_site() {
    let mut log_dis = Array2::from_elem((2, 4), f64::NEG_INFINITY);
    log_dis[[1, 0]] = 0.0;
    log_dis[[1, 1]] = -1.0;
    log_dis[[1, 2]] = -2.0;
    log_dis[[1, 3]] = -3.0;

    let (dis, site_log_norms) = normalize_sites_from_log(log_dis);

    assert_valid_rows(&dis);
    pretty_assert_abs_diff_eq!(array![0.25, 0.25, 0.25, 0.25], dis.row(0), epsilon = 1e-10);
    let sum = 1.0 + (-1.0_f64).exp() + (-2.0_f64).exp() + (-3.0_f64).exp();
    let expected_row1 = array![
      1.0 / sum,
      (-1.0_f64).exp() / sum,
      (-2.0_f64).exp() / sum,
      (-3.0_f64).exp() / sum
    ];
    pretty_assert_abs_diff_eq!(expected_row1, dis.row(1), epsilon = 1e-10);
    pretty_assert_neg_inf!(site_log_norms[0]);
    pretty_assert_abs_diff_eq!(sum.ln(), site_log_norms[1], epsilon = 1e-10);
  }

  #[test]
  fn test_normalize_inplace_divides_rows_by_their_sums_and_returns_sum_of_log_sums() {
    let mut dis = array![[1.0, 2.0, 3.0, 4.0], [4.0, 3.0, 2.0, 1.0]];
    let log_lh = normalize_inplace(&mut dis);

    assert_valid_rows(&dis);
    pretty_assert_abs_diff_eq!(array![[0.1, 0.2, 0.3, 0.4], [0.4, 0.3, 0.2, 0.1]], dis, epsilon = 1e-10);
    pretty_assert_abs_diff_eq!(10.0_f64.ln() + 10.0_f64.ln(), log_lh, epsilon = 1e-10);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::zero(     array![[0.0, 0.0, 0.0, 0.0]])]
  #[case::nan(      array![[f64::NAN, f64::NAN, f64::NAN, f64::NAN]])]
  #[case::infinite( array![[f64::INFINITY, 0.0, 0.0, 0.0]])]
  #[trace]
  fn test_normalize_inplace_replaces_row_without_finite_positive_sum_by_uniform(#[case] dis: Array2<f64>) {
    let mut dis = dis;
    let log_lh = normalize_inplace(&mut dis);

    assert_valid_rows(&dis);
    pretty_assert_abs_diff_eq!(array![[0.25, 0.25, 0.25, 0.25]], dis, epsilon = 1e-10);
    pretty_assert_neg_inf!(log_lh, "log_lh should be NEG_INFINITY for a degenerate row, got {log_lh}");
  }

  #[test]
  fn test_normalize_inplace_returns_neg_inf_when_any_row_is_degenerate() {
    let mut dis = array![[0.0, 0.0, 0.0, 0.0], [2.0, 4.0, 2.0, 2.0]];
    let log_lh = normalize_inplace(&mut dis);

    assert_valid_rows(&dis);
    pretty_assert_abs_diff_eq!(
      array![[0.25, 0.25, 0.25, 0.25], [0.2, 0.4, 0.2, 0.2]],
      dis,
      epsilon = 1e-10
    );
    pretty_assert_neg_inf!(
      log_lh,
      "log_lh should be NEG_INFINITY when any row is zero, got {log_lh}"
    );
  }

  mod helpers {
    use ndarray::Array2;
    use treetime_utils::pretty_assert_abs_diff_eq;

    pub(super) fn assert_valid_rows(dis: &Array2<f64>) {
      for row in dis.rows() {
        assert!(
          row.iter().all(|&v| v >= 0.0 && v.is_finite()),
          "all probabilities must be non-negative and finite: {row}"
        );
        pretty_assert_abs_diff_eq!(1.0, row.sum(), epsilon = 1e-10);
      }
    }
  }
}
