#[cfg(test)]
mod tests {
  use crate::commands::optimize::optimize_unified::{evaluate_mixed, grid_search_branch_lengths};
  use ndarray::array;

  use super::super::test_grid_search_support::tests::{grid_search, make_dense_contribution};

  #[test]
  fn test_grid_search_finds_maximum_log_lh() {
    // Grid search should find the branch length that maximizes log-LH
    let coefficients = array![[0.9, 0.03, 0.03, 0.04], [0.03, 0.9, 0.03, 0.04], [0.1, 0.1, 0.7, 0.1],];
    let contribution = make_dense_contribution(coefficients);
    let contributions = vec![contribution];

    let branch_length = 0.1;
    let one_mutation = 0.001;

    let best_bl = grid_search(&contributions, branch_length, one_mutation);

    // Verify no other point in the grid has higher log-LH
    let branch_lengths = grid_search_branch_lengths(branch_length, one_mutation).unwrap();
    let best_log_lh = evaluate_mixed(&contributions, best_bl).log_lh;

    for &bl in &branch_lengths {
      let log_lh = evaluate_mixed(&contributions, bl).log_lh;
      assert!(
        log_lh <= best_log_lh + 1e-10,
        "Found bl={bl} with log_lh={log_lh} > best_log_lh={best_log_lh}"
      );
    }
  }

  #[test]
  fn test_grid_search_result_in_range() {
    // Result should be within the search range
    let coefficients = array![[0.9, 0.03, 0.03, 0.04], [0.03, 0.9, 0.03, 0.04],];
    let contribution = make_dense_contribution(coefficients);
    let contributions = vec![contribution];

    let branch_length = 0.1;
    let one_mutation = 0.001;

    let best_bl = grid_search(&contributions, branch_length, one_mutation);

    let grid = grid_search_branch_lengths(branch_length, one_mutation).unwrap();
    let lower = grid[0];
    let upper = grid[grid.len() - 1];

    assert!(best_bl >= lower, "best_bl={best_bl} < lower={lower}");
    assert!(best_bl <= upper, "best_bl={best_bl} > upper={upper}");
  }

  #[test]
  fn test_grid_search_with_small_branch_length() {
    // Grid search should work even with very small branch lengths
    let coefficients = array![[0.99, 0.003, 0.003, 0.004], [0.003, 0.99, 0.003, 0.004],];
    let contribution = make_dense_contribution(coefficients);
    let contributions = vec![contribution];

    let branch_length = 0.001;
    let one_mutation = 0.0001;

    let best_bl = grid_search(&contributions, branch_length, one_mutation);

    assert!(best_bl >= 0.0);
    assert!(best_bl.is_finite());
  }

  #[test]
  fn test_grid_search_with_large_branch_length() {
    // Grid search should work with larger branch lengths
    let coefficients = array![[0.5, 0.15, 0.15, 0.2], [0.15, 0.5, 0.15, 0.2],];
    let contribution = make_dense_contribution(coefficients);
    let contributions = vec![contribution];

    let branch_length = 1.0;
    let one_mutation = 0.001;

    let best_bl = grid_search(&contributions, branch_length, one_mutation);

    let grid = grid_search_branch_lengths(branch_length, one_mutation).unwrap();
    let upper = grid[grid.len() - 1];

    assert!(best_bl >= 0.0);
    assert!(best_bl.is_finite());
    assert!(best_bl <= upper, "best_bl={best_bl} > upper={upper}");
  }
}
