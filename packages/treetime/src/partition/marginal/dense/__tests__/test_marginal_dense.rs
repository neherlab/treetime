#[cfg(test)]
mod tests {
  use crate::partition::marginal::shared::{normalize_from_log, normalize_inplace};
  use crate::pretty_assert_neg_inf;
  use approx::assert_abs_diff_eq;
  use ndarray::{Array2, array};

  fn assert_valid_rows(dis: &Array2<f64>) {
    for row in dis.rows() {
      assert!(
        row.iter().all(|&v| v >= 0.0 && v.is_finite()),
        "all probabilities must be non-negative and finite: {row}"
      );
      assert_abs_diff_eq!(row.sum(), 1.0, epsilon = 1e-10);
    }
  }

  #[test]
  fn test_normalize_from_log_equal_probs() {
    let log_dis = Array2::from_elem((2, 4), 0.25_f64.ln());
    let (dis, log_lh) = normalize_from_log(&log_dis);

    assert_valid_rows(&dis);
    assert_abs_diff_eq!(dis, Array2::from_elem((2, 4), 0.25), epsilon = 1e-10);
    assert!(log_lh.is_finite());
  }

  #[test]
  fn test_normalize_from_log_descending() {
    let log_dis = array![[0.0, -1.0, -2.0, -3.0]];
    let (dis, log_lh) = normalize_from_log(&log_dis);

    let sum = 1.0 + (-1.0_f64).exp() + (-2.0_f64).exp() + (-3.0_f64).exp();
    let expected = array![[
      1.0 / sum,
      (-1.0_f64).exp() / sum,
      (-2.0_f64).exp() / sum,
      (-3.0_f64).exp() / sum
    ]];

    assert_valid_rows(&dis);
    assert_abs_diff_eq!(dis, expected, epsilon = 1e-10);
    assert_abs_diff_eq!(log_lh, sum.ln(), epsilon = 1e-10);
  }

  #[test]
  fn test_normalize_from_log_all_neg_inf_single_row() {
    let log_dis = Array2::from_elem((1, 4), f64::NEG_INFINITY);
    let (dis, log_lh) = normalize_from_log(&log_dis);

    assert_valid_rows(&dis);
    assert_abs_diff_eq!(dis, Array2::from_elem((1, 4), 0.25), epsilon = 1e-10);
    pretty_assert_neg_inf!(log_lh, "log_lh should be NEG_INFINITY, got {log_lh}");
  }

  #[test]
  fn test_normalize_from_log_all_neg_inf_multiple_rows() {
    let log_dis = Array2::from_elem((3, 4), f64::NEG_INFINITY);
    let (dis, log_lh) = normalize_from_log(&log_dis);

    assert_valid_rows(&dis);
    assert_abs_diff_eq!(dis, Array2::from_elem((3, 4), 0.25), epsilon = 1e-10);
    pretty_assert_neg_inf!(log_lh, "log_lh should be NEG_INFINITY, got {log_lh}");
  }

  #[test]
  fn test_normalize_from_log_mixed_neg_inf_and_finite_rows() {
    let mut log_dis = Array2::from_elem((2, 4), f64::NEG_INFINITY);
    log_dis[[1, 0]] = 0.0;
    log_dis[[1, 1]] = -1.0;
    log_dis[[1, 2]] = -2.0;
    log_dis[[1, 3]] = -3.0;

    let (dis, log_lh) = normalize_from_log(&log_dis);

    assert_valid_rows(&dis);
    assert_abs_diff_eq!(dis.row(0), array![0.25, 0.25, 0.25, 0.25], epsilon = 1e-10);
    let sum = 1.0 + (-1.0_f64).exp() + (-2.0_f64).exp() + (-3.0_f64).exp();
    let expected_row1 = array![
      1.0 / sum,
      (-1.0_f64).exp() / sum,
      (-2.0_f64).exp() / sum,
      (-3.0_f64).exp() / sum
    ];
    assert_abs_diff_eq!(dis.row(1), expected_row1, epsilon = 1e-10);
    pretty_assert_neg_inf!(log_lh, "log_lh should be NEG_INFINITY, got {log_lh}");
  }

  #[test]
  fn test_normalize_from_log_mixed_finite_neg_inf_within_row() {
    let log_dis = array![[0.0, f64::NEG_INFINITY, -1.0, f64::NEG_INFINITY]];
    let (dis, log_lh) = normalize_from_log(&log_dis);

    let sum = 1.0 + (-1.0_f64).exp();
    let expected = array![[1.0 / sum, 0.0, (-1.0_f64).exp() / sum, 0.0]];

    assert_valid_rows(&dis);
    assert_abs_diff_eq!(dis, expected, epsilon = 1e-10);
    assert_abs_diff_eq!(log_lh, sum.ln(), epsilon = 1e-10);
  }

  #[test]
  fn test_normalize_from_log_large_negative_values() {
    let log_dis = array![[-1000.0, -1001.0, -1002.0, -1003.0]];
    let (dis, log_lh) = normalize_from_log(&log_dis);

    let reference = array![[0.0, -1.0, -2.0, -3.0]];
    let (ref_dis, ref_log_lh) = normalize_from_log(&reference);

    assert_valid_rows(&dis);
    assert_abs_diff_eq!(dis, ref_dis, epsilon = 1e-10);
    assert_abs_diff_eq!(log_lh, ref_log_lh - 1000.0, epsilon = 1e-10);
  }

  #[test]
  fn test_normalize_from_log_three_states() {
    let log_dis = Array2::from_elem((1, 3), f64::NEG_INFINITY);
    let (dis, log_lh) = normalize_from_log(&log_dis);

    assert_valid_rows(&dis);
    assert_abs_diff_eq!(dis, Array2::from_elem((1, 3), 1.0 / 3.0), epsilon = 1e-10);
    pretty_assert_neg_inf!(log_lh, "log_lh should be NEG_INFINITY, got {log_lh}");
  }

  #[test]
  fn test_normalize_inplace_normal_rows() {
    let mut dis = array![[1.0, 2.0, 3.0, 4.0], [4.0, 3.0, 2.0, 1.0]];
    let log_lh = normalize_inplace(&mut dis);

    assert_valid_rows(&dis);
    assert_abs_diff_eq!(dis.row(0), array![0.1, 0.2, 0.3, 0.4], epsilon = 1e-10);
    assert_abs_diff_eq!(dis.row(1), array![0.4, 0.3, 0.2, 0.1], epsilon = 1e-10);
    assert_abs_diff_eq!(log_lh, 10.0_f64.ln() + 10.0_f64.ln(), epsilon = 1e-10);
  }

  #[test]
  fn test_normalize_inplace_zero_row_returns_uniform() {
    let mut dis = array![[0.0, 0.0, 0.0, 0.0]];
    let log_lh = normalize_inplace(&mut dis);

    assert_valid_rows(&dis);
    assert_abs_diff_eq!(dis.row(0), array![0.25, 0.25, 0.25, 0.25], epsilon = 1e-10);
    pretty_assert_neg_inf!(log_lh, "log_lh should be NEG_INFINITY for all-zero row, got {log_lh}");
  }

  #[test]
  fn test_normalize_inplace_mixed_zero_and_normal_rows() {
    let mut dis = array![[0.0, 0.0, 0.0, 0.0], [2.0, 4.0, 2.0, 2.0]];
    let log_lh = normalize_inplace(&mut dis);

    assert_valid_rows(&dis);
    assert_abs_diff_eq!(dis.row(0), array![0.25, 0.25, 0.25, 0.25], epsilon = 1e-10);
    assert_abs_diff_eq!(dis.row(1), array![0.2, 0.4, 0.2, 0.2], epsilon = 1e-10);
    pretty_assert_neg_inf!(
      log_lh,
      "log_lh should be NEG_INFINITY when any row is zero, got {log_lh}"
    );
  }

  #[test]
  fn test_normalize_inplace_nan_row_returns_uniform() {
    let mut dis = array![[f64::NAN, f64::NAN, f64::NAN, f64::NAN]];
    let log_lh = normalize_inplace(&mut dis);

    assert_valid_rows(&dis);
    assert_abs_diff_eq!(dis.row(0), array![0.25, 0.25, 0.25, 0.25], epsilon = 1e-10);
    pretty_assert_neg_inf!(log_lh, "log_lh should be NEG_INFINITY for NaN row, got {log_lh}");
  }

  #[test]
  fn test_normalize_inplace_inf_row_returns_uniform() {
    let mut dis = array![[f64::INFINITY, 0.0, 0.0, 0.0]];
    let log_lh = normalize_inplace(&mut dis);

    assert_valid_rows(&dis);
    assert_abs_diff_eq!(dis.row(0), array![0.25, 0.25, 0.25, 0.25], epsilon = 1e-10);
    pretty_assert_neg_inf!(log_lh, "log_lh should be NEG_INFINITY for inf row, got {log_lh}");
  }
}
