use ndarray::{Array1, Array2, Axis};
use treetime_utils::pretty_assert_abs_diff_eq;

pub(super) fn assert_stochastic_matrix(p: &Array2<f64>, context: &str) {
  let col_sums = p.sum_axis(Axis(0));
  pretty_assert_abs_diff_eq!(col_sums, Array1::ones(p.ncols()), epsilon = 1e-10);
  assert!(
    !p.iter().any(|&x| x < -1e-14),
    "{context}: matrix contains negative value"
  );
}
