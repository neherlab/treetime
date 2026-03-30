#[cfg(test)]
mod tests {
  use ndarray::array;

  use super::super::test_grid_search_support::tests::{grid_search, make_dense_contribution};

  #[test]
  fn test_grid_search_one_mutation_dominates_range() {
    // When one_mutation is large relative to branch_length,
    // the range should still be valid
    let coefficients = array![[0.9, 0.03, 0.03, 0.04], [0.03, 0.9, 0.03, 0.04],];
    let contribution = make_dense_contribution(coefficients);
    let contributions = vec![contribution];

    let branch_length = 0.001;
    let one_mutation = 0.01; // 10x branch_length

    let upper = 1.5 * branch_length + one_mutation;

    let best_bl = grid_search(&contributions, branch_length, one_mutation);
    assert!(best_bl >= 0.0);
    assert!(best_bl <= upper);
  }

  #[test]
  fn test_grid_search_returns_non_negative_branch_length() {
    // Grid search returns non-negative branch length (zero is valid when
    // zero is the likelihood maximum, e.g. zero-mutation edges)
    let coefficients = array![[0.9, 0.03, 0.03, 0.04], [0.03, 0.9, 0.03, 0.04],];
    let contribution = make_dense_contribution(coefficients);
    let contributions = vec![contribution];

    let test_cases = [(0.1, 0.001), (0.01, 0.001), (0.001, 0.0001), (1.0, 0.001)];

    for (branch_length, one_mutation) in test_cases {
      let best_bl = grid_search(&contributions, branch_length, one_mutation);
      assert!(
        best_bl >= 0.0,
        "grid_search({branch_length}, {one_mutation}) = {best_bl} should be >= 0"
      );
    }
  }
}
