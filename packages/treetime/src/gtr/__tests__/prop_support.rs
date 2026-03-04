use ndarray::{Array1, Array2};
use proptest::test_runner::TestCaseError;

/// Assert each column of `matrix` sums to `expected` within `epsilon`.
pub fn prop_assert_columns_sum_to(
  matrix: &Array2<f64>,
  expected: f64,
  epsilon: f64,
) -> Result<(), TestCaseError> {
  for (j, col) in matrix.columns().into_iter().enumerate() {
    let sum = col.sum();
    if !approx::abs_diff_eq!(sum, expected, epsilon = epsilon) {
      return Err(TestCaseError::fail(format!(
        "column {j}: sum = {sum}, expected {expected} (epsilon = {epsilon})"
      )));
    }
  }
  Ok(())
}

/// Assert each row of `matrix` sums to `expected` within `epsilon`.
pub fn prop_assert_rows_sum_to(
  matrix: &Array2<f64>,
  expected: f64,
  epsilon: f64,
) -> Result<(), TestCaseError> {
  for (i, row) in matrix.rows().into_iter().enumerate() {
    let sum = row.sum();
    if !approx::abs_diff_eq!(sum, expected, epsilon = epsilon) {
      return Err(TestCaseError::fail(format!(
        "row {i}: sum = {sum}, expected {expected} (epsilon = {epsilon})"
      )));
    }
  }
  Ok(())
}

/// Assert detailed balance: pi[j] * Q[i,j] = pi[i] * Q[j,i] for all i != j.
pub fn prop_assert_detailed_balance(
  q: &Array2<f64>,
  pi: &Array1<f64>,
  epsilon: f64,
) -> Result<(), TestCaseError> {
  let n = q.nrows();
  for i in 0..n {
    for j in 0..n {
      if i != j {
        let flux_ji = pi[j] * q[[i, j]];
        let flux_ij = pi[i] * q[[j, i]];
        if !approx::abs_diff_eq!(flux_ji, flux_ij, epsilon = epsilon) {
          return Err(TestCaseError::fail(format!(
            "detailed balance violated: pi[{j}]*Q[{i},{j}] = {flux_ji}, pi[{i}]*Q[{j},{i}] = {flux_ij} (epsilon = {epsilon})"
          )));
        }
      }
    }
  }
  Ok(())
}
