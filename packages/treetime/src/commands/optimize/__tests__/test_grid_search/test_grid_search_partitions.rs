use crate::commands::optimize::optimize_unified::evaluate_mixed;
use ndarray::{Array1, array};

use super::test_grid_search_support::{grid_search, make_dense_contribution};

#[test]
fn test_grid_search_multiple_partitions() {
  // Grid search should work with multiple contributions
  let coefficients1 = array![[0.9, 0.03, 0.03, 0.04], [0.03, 0.9, 0.03, 0.04],];
  let coefficients2 = array![[0.8, 0.06, 0.06, 0.08], [0.06, 0.8, 0.06, 0.08],];
  let contribution1 = make_dense_contribution(coefficients1);
  let contribution2 = make_dense_contribution(coefficients2);
  let contributions = vec![contribution1, contribution2];

  let branch_length = 0.1;
  let one_mutation = 0.001;

  let best_bl = grid_search(&contributions, branch_length, one_mutation);

  // Verify combined log-LH is maximized
  let branch_lengths = Array1::linspace(0.1 * one_mutation, 1.5 * branch_length + one_mutation, 100);
  let best_log_lh = evaluate_mixed(&contributions, best_bl).log_lh;

  for &bl in &branch_lengths {
    let log_lh = evaluate_mixed(&contributions, bl).log_lh;
    assert!(log_lh <= best_log_lh + 1e-10);
  }
}

#[test]
fn test_grid_search_uniform_coefficients() {
  // With uniform coefficients, any branch length gives similar log-LH
  // Grid search should still return a valid result
  let coefficients = array![[0.25, 0.25, 0.25, 0.25], [0.25, 0.25, 0.25, 0.25],];
  let contribution = make_dense_contribution(coefficients);
  let contributions = vec![contribution];

  let branch_length = 0.1;
  let one_mutation = 0.001;

  let best_bl = grid_search(&contributions, branch_length, one_mutation);

  assert!(best_bl.is_finite());
  assert!(best_bl >= 0.0);
}
