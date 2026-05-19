use ndarray::{Array1, Array2, Axis};
use proptest::test_runner::TestCaseError;
use treetime_utils::prop_assert_array_abs_diff_eq;

/// Assert each column of `matrix` sums to `expected` within `epsilon`.
pub fn prop_assert_columns_sum_to(matrix: &Array2<f64>, expected: f64, epsilon: f64) -> Result<(), TestCaseError> {
  let col_sums = matrix.sum_axis(Axis(0));
  let expected_arr = Array1::from_elem(matrix.ncols(), expected);
  prop_assert_array_abs_diff_eq!(col_sums, expected_arr, epsilon = epsilon);
  Ok(())
}

/// Assert each row of `matrix` sums to `expected` within `epsilon`.
pub fn prop_assert_rows_sum_to(matrix: &Array2<f64>, expected: f64, epsilon: f64) -> Result<(), TestCaseError> {
  let row_sums = matrix.sum_axis(Axis(1));
  let expected_arr = Array1::from_elem(matrix.nrows(), expected);
  prop_assert_array_abs_diff_eq!(row_sums, expected_arr, epsilon = epsilon);
  Ok(())
}

/// Assert detailed balance: pi[j] * Q[i,j] = pi[i] * Q[j,i] for all i != j.
///
/// Equivalent to checking that the flux matrix F[i,j] = pi[j] * Q[i,j] is symmetric.
pub fn prop_assert_detailed_balance(q: &Array2<f64>, pi: &Array1<f64>, epsilon: f64) -> Result<(), TestCaseError> {
  let flux = q * &pi.view().insert_axis(Axis(0));
  let flux_t = flux.t().to_owned();
  prop_assert_array_abs_diff_eq!(flux, flux_t, epsilon = epsilon);
  Ok(())
}
