#[cfg(test)]
mod tests {
  use crate::commands::optimize::optimize_unified::grid_search_branch_lengths;
  use approx::assert_abs_diff_eq;
  use ndarray::Array2;

  use super::super::test_grid_search_support::tests::{grid_search, make_dense_contribution};

  #[test]
  fn test_grid_search_branch_lengths_zero_branch_covers_full_range() {
    // The original bug: when branch_length=0.0 and L=1000, the grid covered
    // only [1e-4, 1e-3]. After the fix, the grid must extend to at least 0.5
    // subs/site regardless of the current branch length.
    let one_mutation = 1.0 / 1000.0; // L = 1000
    let grid = grid_search_branch_lengths(0.0, one_mutation);

    let lower = grid[0];
    let upper = grid[grid.len() - 1];

    // Lower bound: 0.1 * one_mutation = 1e-4
    assert_abs_diff_eq!(lower, 0.1 * one_mutation, epsilon = 1e-10);
    // Upper bound: max(1.5*0.0 + 1e-3, 0.5) = 0.5
    assert_abs_diff_eq!(upper, 0.5, epsilon = 1e-10);
    // Grid must span at least 3 orders of magnitude
    assert!(
      upper / lower > 1000.0,
      "Grid range {lower}..{upper} spans less than 3 orders of magnitude"
    );
  }

  #[test]
  fn test_grid_search_branch_lengths_large_branch_extends_beyond_minimum() {
    // When current branch is longer than GRID_SEARCH_MIN_UPPER, the proportional
    // bound (1.5 * branch_length + one_mutation) should dominate.
    let one_mutation = 0.001;
    let branch_length = 1.0;
    let grid = grid_search_branch_lengths(branch_length, one_mutation);

    let upper = grid[grid.len() - 1];
    let expected_upper = 1.5 * branch_length + one_mutation;

    assert_abs_diff_eq!(upper, expected_upper, epsilon = 1e-10);
  }

  #[test]
  fn test_grid_search_branch_lengths_is_log_spaced() {
    // Log-spaced grid has monotonically increasing spacing: the gap between
    // consecutive points grows as the branch length increases.
    let grid = grid_search_branch_lengths(0.0, 0.001);

    assert!(grid.len() == 100);

    // All values must be positive and strictly increasing
    for i in 1..grid.len() {
      assert!(
        grid[i] > grid[i - 1],
        "Grid not strictly increasing at index {i}: {} >= {}",
        grid[i - 1],
        grid[i]
      );
    }

    // Spacing must increase (log-spacing property)
    let first_gap = grid[1] - grid[0];
    let last_gap = grid[grid.len() - 1] - grid[grid.len() - 2];
    assert!(
      last_gap > first_gap * 10.0,
      "Last gap ({last_gap}) should be much larger than first gap ({first_gap}) for log-spaced grid"
    );

    // Ratios between consecutive points should be approximately constant
    let ratio_start = grid[1] / grid[0];
    let ratio_end = grid[grid.len() - 1] / grid[grid.len() - 2];
    assert_abs_diff_eq!(ratio_start, ratio_end, epsilon = 1e-10);
  }

  #[test]
  fn test_grid_search_branch_lengths_has_correct_point_count() {
    let grid = grid_search_branch_lengths(0.0, 0.001);
    assert_eq!(grid.len(), 100);

    let grid = grid_search_branch_lengths(0.5, 0.0001);
    assert_eq!(grid.len(), 100);
  }

  #[test]
  fn test_grid_search_reaches_extended_range_from_zero() {
    // JC69 eigh(Lower) returns eigenvalues in ascending order:
    //   eigvals = [-4/3, -4/3, -4/3, 0]
    // The zero eigenvalue (stationary component) is at index 3.
    // Per-site likelihood: L_i(t) = (k0+k1+k2)*exp(lambda*t) + k3
    //
    // Negative sum (k0+k1+k2 < 0) indicates substitution evidence:
    // L_i(t) increases monotonically with t. For all-substitution
    // sites, the grid search returns the highest grid point.
    //
    // With the old narrow grid at branch_length=0 (L=1000), the upper
    // bound was ~0.001. The fix extends it to 0.5.
    let coefficients = Array2::from_shape_fn((3, 4), |(_, j)| {
      // Negative non-stationary (indices 0-2), positive stationary (index 3)
      [-0.06, -0.06, -0.06, 0.25][j]
    });
    let contribution = make_dense_contribution(coefficients);
    let contributions = vec![contribution];

    let branch_length = 0.0;
    let one_mutation = 0.001; // L = 1000

    let best_bl = grid_search(&contributions, branch_length, one_mutation);

    // With monotonically increasing likelihood, grid search returns near
    // the upper bound (~0.5). The old grid would have returned ~0.001.
    assert!(
      best_bl > 0.1,
      "Expected best_bl > 0.1 (grid should extend to 0.5), got {best_bl}"
    );
    assert!(best_bl.is_finite());
  }

  #[test]
  fn test_grid_search_branch_lengths_coverage_invariant() {
    // For any branch_length >= 0 and any reasonable one_mutation, the grid
    // must always cover at least [one_mutation, 0.5].
    #[rustfmt::skip]
    let cases: &[(f64, f64)] = &[
      (0.0,    0.001  ),
      (0.0,    0.0001 ),
      (0.0,    0.01   ),
      (0.001,  0.001  ),
      (0.1,    0.001  ),
      (0.5,    0.001  ),
      (1.0,    0.001  ),
      (2.0,    0.0001 ),
    ];

    for &(branch_length, one_mutation) in cases {
      let grid = grid_search_branch_lengths(branch_length, one_mutation);
      let lower = grid[0];
      let upper = grid[grid.len() - 1];

      assert!(
        lower <= one_mutation,
        "branch_length={branch_length}, one_mutation={one_mutation}: lower={lower} > one_mutation"
      );
      assert!(
        upper >= 0.5_f64.min(1.5 * branch_length + one_mutation),
        "branch_length={branch_length}, one_mutation={one_mutation}: upper={upper} < min(0.5, proportional)"
      );
    }
  }
}
