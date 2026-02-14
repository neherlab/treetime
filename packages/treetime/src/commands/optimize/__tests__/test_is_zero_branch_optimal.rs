#[cfg(test)]
mod tests {
  use crate::commands::optimize::optimize_dense;
  use crate::commands::optimize::optimize_sparse;
  use crate::commands::optimize::optimize_unified::{OptimizationContribution, is_zero_branch_optimal};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use ndarray::{Array2, array};

  /// Create a dense contribution with specified coefficients and JC69 GTR model.
  ///
  /// For JC69 on nucleotides:
  /// - eigvals = [0, -4/3, -4/3, -4/3] (first eigenvalue is always 0)
  /// - At branch_length=0, likelihood = product of coefficient row sums
  /// - Derivative at zero = sum over positions of (eigval-weighted coeff sum / coeff sum)
  fn make_dense_contribution(coefficients: Array2<f64>) -> OptimizationContribution {
    let gtr = jc69(JC69Params::default()).unwrap();
    OptimizationContribution::Dense(optimize_dense::PartitionContribution::new(coefficients, gtr))
  }

  /// Create a sparse contribution with specified site contributions and JC69 eigenvalues.
  fn make_sparse_contribution(sites: Vec<(f64, Vec<f64>)>) -> OptimizationContribution {
    let gtr = jc69(JC69Params::default()).unwrap();
    let eigenvalues = gtr.eigvals;

    let site_contributions = sites
      .into_iter()
      .map(|(multiplicity, coeffs)| optimize_sparse::SiteContribution {
        multiplicity,
        coefficients: ndarray::Array1::from_vec(coeffs),
      })
      .collect();

    OptimizationContribution::Sparse(optimize_sparse::PartitionContribution {
      site_contributions,
      eigenvalues,
    })
  }

  #[test]
  fn test_is_zero_branch_optimal_empty_contributions() {
    let contributions: Vec<OptimizationContribution> = vec![];

    let result = is_zero_branch_optimal(&contributions);

    // Empty product is 1.0 (>0.01), empty sum is 0.0 (not <0.0), so false
    assert!(!result);
  }

  #[test]
  fn test_is_zero_branch_optimal_high_lh_negative_derivative_returns_true() {
    // Coefficients where:
    // - Zero LH is high (row sums product > 0.01)
    // - Derivative at zero is negative
    //
    // For a single position with coefficients [1, 0, 0, 0]:
    // - LH at zero = sum = 1.0
    // - Derivative = (eigval[0]*1 + eigval[1]*0 + ...) / 1.0 = 0 (since eigval[0]=0)
    //
    // To get negative derivative, we need non-zero coefficients on negative eigenvalues.
    // JC69 eigvals: [0, -4/3, -4/3, -4/3]
    // Coefficients [0, 1, 0, 0] gives derivative = -4/3 / 1 = -4/3 < 0
    let coefficients = array![[0.0, 1.0, 0.0, 0.0]];
    let contribution = make_dense_contribution(coefficients);

    let result = is_zero_branch_optimal(&[contribution]);

    assert!(result);
  }

  #[test]
  fn test_is_zero_branch_optimal_high_lh_positive_derivative_returns_false() {
    // Coefficients where derivative is non-negative.
    // JC69 eigenvalues are sorted ascending, with the zero eigenvalue last.
    // Putting weight only on the last eigenvalue (index 3, closest to 0) gives derivative ~0.
    // For JC69: eigvals ~ [-4/3, -4/3, -4/3, 0] (ascending order)
    let coefficients = array![[0.0, 0.0, 0.0, 1.0]];
    let contribution = make_dense_contribution(coefficients);

    let result = is_zero_branch_optimal(&[contribution]);

    assert!(!result);
  }

  #[test]
  fn test_is_zero_branch_optimal_low_lh_skips_derivative() {
    // Even if derivative would be negative, low LH means we skip the check.
    // Very small coefficients give LH < 0.01
    let coefficients = array![[0.0, 0.001, 0.0, 0.0]];
    let contribution = make_dense_contribution(coefficients);

    // LH = 0.001 < 0.01, so derivative check is skipped
    let result = is_zero_branch_optimal(&[contribution]);

    assert!(!result);
  }

  #[test]
  fn test_is_zero_branch_optimal_multiple_partitions_product_lh() {
    // Two contributions with LH=0.5 each, combined LH = 0.25 > 0.01
    // Both with negative derivatives
    let coeff1 = array![[0.0, 0.5, 0.0, 0.0]];
    let coeff2 = array![[0.0, 0.5, 0.0, 0.0]];
    let contrib1 = make_dense_contribution(coeff1);
    let contrib2 = make_dense_contribution(coeff2);

    let result = is_zero_branch_optimal(&[contrib1, contrib2]);

    // LH = 0.5 * 0.5 = 0.25 > 0.01
    // Derivative = -4/3 + -4/3 = -8/3 < 0
    assert!(result);
  }

  #[test]
  fn test_is_zero_branch_optimal_multiple_partitions_low_combined_lh() {
    // Two contributions that individually have ok LH but product is low
    let coeff1 = array![[0.0, 0.05, 0.0, 0.0]];
    let coeff2 = array![[0.0, 0.05, 0.0, 0.0]];
    let contrib1 = make_dense_contribution(coeff1);
    let contrib2 = make_dense_contribution(coeff2);

    let result = is_zero_branch_optimal(&[contrib1, contrib2]);

    // LH = 0.05 * 0.05 = 0.0025 < 0.01, so false
    assert!(!result);
  }

  #[test]
  fn test_is_zero_branch_optimal_sparse_high_lh_negative_derivative() {
    // Sparse contribution with single variable site
    // multiplicity=1.0, coefficients with weight on negative eigenvalue
    let contribution = make_sparse_contribution(vec![(1.0, vec![0.0, 1.0, 0.0, 0.0])]);

    let result = is_zero_branch_optimal(&[contribution]);

    // LH = 1.0^1 = 1.0 > 0.01
    // Derivative at zero uses weighted coefficients over eigenvalues
    assert!(result);
  }

  #[test]
  fn test_is_zero_branch_optimal_sparse_high_lh_zero_derivative() {
    // Sparse contribution where derivative is zero (all weight on eigval=0)
    // JC69 eigenvalues sorted ascending, zero eigenvalue is last (index 3)
    let contribution = make_sparse_contribution(vec![(1.0, vec![0.0, 0.0, 0.0, 1.0])]);

    let result = is_zero_branch_optimal(&[contribution]);

    // LH = 1.0 > 0.01
    // Derivative = 0 (not < 0)
    assert!(!result);
  }

  #[test]
  fn test_is_zero_branch_optimal_mixed_dense_and_sparse() {
    // Mix dense and sparse contributions
    let dense = make_dense_contribution(array![[0.0, 0.5, 0.0, 0.0]]);
    let sparse = make_sparse_contribution(vec![(1.0, vec![0.0, 0.5, 0.0, 0.0])]);

    let result = is_zero_branch_optimal(&[dense, sparse]);

    // LH = 0.5 * 0.5 = 0.25 > 0.01
    // Both have negative derivatives
    assert!(result);
  }

  #[test]
  fn test_is_zero_branch_optimal_multiple_positions_dense() {
    // Multiple positions in a single dense contribution
    // Each row is one position
    let coefficients = array![
      [0.0, 1.0, 0.0, 0.0], // position 1: negative derivative
      [0.0, 1.0, 0.0, 0.0], // position 2: negative derivative
    ];
    let contribution = make_dense_contribution(coefficients);

    let result = is_zero_branch_optimal(&[contribution]);

    // LH = 1.0 * 1.0 = 1.0 > 0.01
    // Derivative = -4/3 + -4/3 = -8/3 < 0
    assert!(result);
  }

  #[test]
  fn test_is_zero_branch_optimal_sparse_multiplicity_effect() {
    // Multiplicity affects how site contributes to likelihood
    // multiplicity=2 means coeff.sum()^2
    let contribution = make_sparse_contribution(vec![(2.0, vec![0.0, 0.5, 0.0, 0.0])]);

    let result = is_zero_branch_optimal(&[contribution]);

    // LH = 0.5^2 = 0.25 > 0.01
    // Derivative is also affected by multiplicity
    assert!(result);
  }

  #[test]
  fn test_is_zero_branch_optimal_boundary_lh_just_above_threshold() {
    // LH just above 0.01 threshold
    let coefficients = array![[0.0, 0.011, 0.0, 0.0]];
    let contribution = make_dense_contribution(coefficients);

    let result = is_zero_branch_optimal(&[contribution]);

    // LH = 0.011 > 0.01, derivative < 0
    assert!(result);
  }

  #[test]
  fn test_is_zero_branch_optimal_boundary_lh_at_threshold() {
    // LH exactly at 0.01 threshold (condition is >0.01, not >=)
    let coefficients = array![[0.0, 0.01, 0.0, 0.0]];
    let contribution = make_dense_contribution(coefficients);

    let result = is_zero_branch_optimal(&[contribution]);

    // LH = 0.01, not > 0.01, so false
    assert!(!result);
  }
}
