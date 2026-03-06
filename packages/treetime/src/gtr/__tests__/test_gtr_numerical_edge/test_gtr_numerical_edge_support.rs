//! Numerical stability edge case tests for GTR models.
//!
//! Tests boundary conditions where floating-point issues commonly arise:
//! branch length extremes, kappa extremes, skewed frequencies, mu scaling.
//!
//! NOTE: This implementation uses a transposed convention where P(t) is
//! column-stochastic (columns sum to 1), not row-stochastic.

use approx::assert_abs_diff_eq;
use ndarray::{Array1, Array2, Axis};

/// Assert matrix is column-stochastic: columns sum to 1, entries non-negative.
/// Uses transposed convention where P[i, j] = Pr(to state i | from state j).
pub(super) fn assert_stochastic_matrix(p: &Array2<f64>, context: &str) {
  let col_sums = p.sum_axis(Axis(0));
  assert_abs_diff_eq!(col_sums, Array1::ones(p.ncols()), epsilon = 1e-10);
  assert!(
    !p.iter().any(|&x| x < -1e-14),
    "{context}: matrix contains negative value"
  );
}
