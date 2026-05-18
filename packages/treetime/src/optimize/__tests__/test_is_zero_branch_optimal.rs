#[cfg(test)]
mod tests {
  use crate::gtr::get_gtr::{
    F81Params, HKY85Params, JC69Params, Jtt92Params, K80Params, T92Params, TN93Params, f81, hky85, jc69, jtt92, k80,
    t92, tn93,
  };
  use crate::gtr::gtr::GTR;
  use crate::optimize::likelihood::evaluate_mixed_log_lh_only;
  use crate::optimize::zero_boundary::is_zero_branch_optimal;
  use crate::representation::partition::optimization_contribution::OptimizationContribution;
  use crate::representation::partition::optimize_dense;
  use crate::representation::partition::optimize_sparse;
  use approx::assert_abs_diff_eq;
  use ndarray::{Array2, array};
  use rstest::rstest;

  /// Create a dense contribution with specified coefficients and JC69 GTR model.
  ///
  /// JC69 nucleotide eigenvalues (sorted ascending): [-4/3, -4/3, -4/3, 0].
  /// At t=0, per-site likelihood L_i(0) = row sum of coefficients.
  /// Per-site derivative ratio = (eigval-weighted row sum) / (row sum).
  fn make_dense_contribution(coefficients: Array2<f64>) -> OptimizationContribution {
    let gtr = jc69(JC69Params::default()).unwrap();
    OptimizationContribution::Dense(optimize_dense::PartitionContribution::new(coefficients, gtr))
  }

  /// Create a sparse contribution with specified (multiplicity, coefficients) pairs.
  fn make_sparse_contribution(sites: Vec<(f64, Vec<f64>)>) -> OptimizationContribution {
    let gtr = jc69(JC69Params::default()).unwrap();

    let site_contributions = sites
      .into_iter()
      .map(|(multiplicity, coeffs)| optimize_sparse::SiteContribution {
        multiplicity,
        coefficients: ndarray::Array1::from_vec(coeffs),
      })
      .collect();

    OptimizationContribution::Sparse(optimize_sparse::PartitionContribution {
      site_contributions,
      gtr,
    })
  }

  /// Empty contributions: no sites, derivative sum is 0.0 (not < 0).
  #[test]
  fn test_is_zero_branch_optimal_empty_contributions() {
    assert!(!is_zero_branch_optimal(&[]));
  }

  /// Negative derivative at t=0 means zero is optimal.
  /// Coefficients [0, 1, 0, 0]: L_i(0) = 1, derivative = -4/3 < 0.
  #[test]
  fn test_is_zero_branch_optimal_negative_derivative_returns_true() {
    let contribution = make_dense_contribution(array![[0.0, 1.0, 0.0, 0.0]]);
    assert!(is_zero_branch_optimal(&[contribution]));
  }

  /// Zero derivative at t=0 means zero is not a local maximum.
  /// Coefficients [0, 0, 0, 1]: all weight on the zero eigenvalue.
  /// L_i(0) = 1, derivative = 0 * 1 / 1 = 0 (not < 0).
  #[test]
  fn test_is_zero_branch_optimal_zero_derivative_returns_false() {
    let contribution = make_dense_contribution(array![[0.0, 0.0, 0.0, 1.0]]);
    assert!(!is_zero_branch_optimal(&[contribution]));
  }

  /// The decision must be the same for one site vs many identical sites.
  /// The old raw-product approach would underflow to 0.0 for many sites,
  /// flipping the decision. The derivative-sign approach is additive over
  /// sites, so the sign is preserved.
  #[test]
  fn test_is_zero_branch_optimal_scale_invariant_dense() {
    let one_site = make_dense_contribution(array![[0.0, 1.0, 0.0, 0.0]]);
    let result_one = is_zero_branch_optimal(&[one_site]);

    // 100 identical sites: same derivative sign, just larger magnitude
    let many_rows: Vec<[f64; 4]> = vec![[0.0, 1.0, 0.0, 0.0]; 100];
    let arr = Array2::from(many_rows);
    let many_sites = make_dense_contribution(arr);
    let result_many = is_zero_branch_optimal(&[many_sites]);

    assert_eq!(
      result_one, result_many,
      "decision must not depend on number of identical sites"
    );
  }

  /// Scale invariance for sparse: same single-site pattern repeated via
  /// multiplicity vs duplicated as separate sites.
  #[test]
  fn test_is_zero_branch_optimal_scale_invariant_sparse() {
    let single = make_sparse_contribution(vec![(1.0, vec![0.0, 1.0, 0.0, 0.0])]);
    let result_single = is_zero_branch_optimal(&[single]);

    let high_mult = make_sparse_contribution(vec![(100.0, vec![0.0, 1.0, 0.0, 0.0])]);
    let result_mult = is_zero_branch_optimal(&[high_mult]);

    assert_eq!(
      result_single, result_mult,
      "decision must not depend on multiplicity scaling"
    );
  }

  /// Many ambiguous sites with L_i(0) = 0.25 each. The old implementation
  /// would compute 0.25^N which underflows to 0.0 for large N, preventing
  /// the derivative check. The new implementation checks each site
  /// independently and evaluates the derivative sign.
  ///
  /// JC69 eigvals [-4/3, -4/3, -4/3, 0]. For uniform coefficients
  /// [0.0625, 0.0625, 0.0625, 0.0625], derivative ratio per site =
  /// (-4/3 * 3 * 0.0625) / (4 * 0.0625) = -1.0. Negative, so zero is optimal.
  #[test]
  fn test_is_zero_branch_optimal_underflow_resistant() {
    let many_rows: Vec<[f64; 4]> = vec![[0.0625, 0.0625, 0.0625, 0.0625]; 1000];
    let arr = Array2::from(many_rows);
    let contribution = make_dense_contribution(arr);

    // Old code: 0.25^1000 = 0.0, would return false. New code: derivative = -1000 < 0.
    assert!(is_zero_branch_optimal(&[contribution]));
  }

  /// A site with all-zero coefficients has L_i(0) = 0. The derivative
  /// formula would divide by zero. The shortcut must decline to decide.
  #[test]
  fn test_is_zero_branch_optimal_zero_site_lh_returns_false() {
    let coefficients = array![
      [0.0, 1.0, 0.0, 0.0], // valid site
      [0.0, 0.0, 0.0, 0.0], // degenerate: L_i(0) = 0
    ];
    let contribution = make_dense_contribution(coefficients);
    assert!(!is_zero_branch_optimal(&[contribution]));
  }

  /// A site with negative coefficient sum (physically impossible but
  /// numerically possible from floating-point error) must be rejected.
  #[test]
  fn test_is_zero_branch_optimal_negative_site_lh_returns_false() {
    let coefficients = array![[1.0, -2.0, 0.0, 0.0]];
    let contribution = make_dense_contribution(coefficients);
    assert!(!is_zero_branch_optimal(&[contribution]));
  }

  /// A site with non-finite coefficient sum (NaN or Inf) must be rejected.
  #[test]
  fn test_is_zero_branch_optimal_nonfinite_site_lh_returns_false() {
    let coefficients = array![[f64::INFINITY, 0.0, 0.0, 0.0]];
    let contribution = make_dense_contribution(coefficients);
    assert!(!is_zero_branch_optimal(&[contribution]));
  }

  /// Two contributions both with negative derivatives combine additively.
  #[test]
  fn test_is_zero_branch_optimal_multiple_partitions_negative_derivative() {
    let contrib1 = make_dense_contribution(array![[0.0, 0.5, 0.0, 0.0]]);
    let contrib2 = make_dense_contribution(array![[0.0, 0.5, 0.0, 0.0]]);
    assert!(is_zero_branch_optimal(&[contrib1, contrib2]));
  }

  /// Small per-site likelihoods with negative derivative. The old
  /// raw-product approach would compute 0.05 * 0.05 = 0.0025 < 0.01 and
  /// skip the derivative check. The new approach checks sites individually
  /// (both are positive and finite) and evaluates derivative sign.
  #[test]
  fn test_is_zero_branch_optimal_small_lh_still_decides() {
    let contrib1 = make_dense_contribution(array![[0.0, 0.05, 0.0, 0.0]]);
    let contrib2 = make_dense_contribution(array![[0.0, 0.05, 0.0, 0.0]]);

    // Per-site L_i(0) = 0.05 > 0: valid. Derivative = -4/3 + -4/3 < 0.
    assert!(is_zero_branch_optimal(&[contrib1, contrib2]));
  }

  /// Mixed dense and sparse contributions with negative derivatives.
  #[test]
  fn test_is_zero_branch_optimal_mixed_dense_and_sparse() {
    let dense = make_dense_contribution(array![[0.0, 0.5, 0.0, 0.0]]);
    let sparse = make_sparse_contribution(vec![(1.0, vec![0.0, 0.5, 0.0, 0.0])]);
    assert!(is_zero_branch_optimal(&[dense, sparse]));
  }

  /// Multiple positions in one dense contribution. Derivative sums over
  /// all positions: -4/3 + -4/3 = -8/3 < 0.
  #[test]
  fn test_is_zero_branch_optimal_multiple_positions_dense() {
    #[rustfmt::skip]
    let coefficients = array![
      [0.0, 1.0, 0.0, 0.0],
      [0.0, 1.0, 0.0, 0.0],
    ];
    let contribution = make_dense_contribution(coefficients);
    assert!(is_zero_branch_optimal(&[contribution]));
  }

  /// Sparse site with negative derivative.
  #[test]
  fn test_is_zero_branch_optimal_sparse_negative_derivative() {
    let contribution = make_sparse_contribution(vec![(1.0, vec![0.0, 1.0, 0.0, 0.0])]);
    assert!(is_zero_branch_optimal(&[contribution]));
  }

  /// Sparse site where derivative is zero (all weight on zero eigenvalue).
  #[test]
  fn test_is_zero_branch_optimal_sparse_zero_derivative() {
    let contribution = make_sparse_contribution(vec![(1.0, vec![0.0, 0.0, 0.0, 1.0])]);
    assert!(!is_zero_branch_optimal(&[contribution]));
  }

  /// Sparse multiplicity scales the derivative but not the sign.
  #[test]
  fn test_is_zero_branch_optimal_sparse_multiplicity_preserves_sign() {
    let contribution = make_sparse_contribution(vec![(50.0, vec![0.0, 0.5, 0.0, 0.0])]);
    assert!(is_zero_branch_optimal(&[contribution]));
  }

  /// Coefficients that produce a positive finite L_i(0) but overflow in the
  /// derivative ratio. This exercises the `derivative.is_finite()` guard at
  /// zero_boundary.rs which would otherwise be a dead code path.
  ///
  /// With JC69's equal non-zero eigenvalues (-4/3), the per-site derivative
  /// ratio is always bounded between -4/3 and 0, so overflow is impossible.
  /// To trigger overflow, we need a GTR model with a very large eigenvalue
  /// combined with near-canceling coefficients:
  ///
  /// eigvals = [-1e300, 0, 0, 0], coefficients = [1, 0, 0, -1 + 2*eps]:
  /// - L_i(0) = 2 * f64::EPSILON ≈ 4.44e-16 (positive and finite)
  /// - numerator = 1 * (-1e300) = -1e300
  /// - ratio = -1e300 / 4.44e-16 ≈ -2.25e315 (overflows f64::MAX)
  #[test]
  fn test_is_zero_branch_optimal_nonfinite_derivative_returns_false() {
    let mut gtr = jc69(JC69Params::default()).unwrap();
    gtr.eigvals = array![-1e300, 0.0, 0.0, 0.0];
    let coefficients = array![[1.0, 0.0, 0.0, -1.0 + f64::EPSILON * 2.0]];
    let contribution = OptimizationContribution::Dense(optimize_dense::PartitionContribution::new(coefficients, gtr));

    // Site likelihood is positive and finite (2 * f64::EPSILON)
    assert!(contribution.all_sites_valid_at_zero());
    // But the derivative overflows due to large eigenvalue / tiny denominator
    assert!(!contribution.zero_branch_length_derivative().is_finite());
    // So is_zero_branch_optimal declines to decide
    assert!(!is_zero_branch_optimal(&[contribution]));
  }

  /// Verify the exact numerical value of the derivative, not just its sign.
  ///
  /// JC69 eigvals [-4/3, -4/3, -4/3, 0]. Coefficients [0, 1, 0, 0]:
  /// - L_i(0) = 1
  /// - numerator = 0*(-4/3) + 1*(-4/3) + 0*(-4/3) + 0*0 = -4/3
  /// - derivative = -4/3 / 1 = -4/3
  #[test]
  fn test_is_zero_branch_optimal_derivative_magnitude_dense() {
    let contribution = make_dense_contribution(array![[0.0, 1.0, 0.0, 0.0]]);
    let derivative = contribution.zero_branch_length_derivative();
    assert_abs_diff_eq!(-4.0 / 3.0, derivative, epsilon = 1e-10);
  }

  /// Sparse derivative with multiplicity scales the magnitude.
  /// multiplicity=3, same coefficients: derivative = 3 * (-4/3) = -4.
  #[test]
  fn test_is_zero_branch_optimal_derivative_magnitude_sparse() {
    let contribution = make_sparse_contribution(vec![(3.0, vec![0.0, 1.0, 0.0, 0.0])]);
    let derivative = contribution.zero_branch_length_derivative();
    assert_abs_diff_eq!(-4.0, derivative, epsilon = 1e-10);
  }

  /// Multi-site derivative sums per-site values.
  /// Two sites each contributing -4/3: total = -8/3.
  #[test]
  fn test_is_zero_branch_optimal_derivative_magnitude_multi_site() {
    #[rustfmt::skip]
    let coefficients = array![
      [0.0, 1.0, 0.0, 0.0],
      [0.0, 1.0, 0.0, 0.0],
    ];
    let contribution = make_dense_contribution(coefficients);
    let derivative = contribution.zero_branch_length_derivative();
    assert_abs_diff_eq!(-8.0 / 3.0, derivative, epsilon = 1e-10);
  }

  /// Dense and sparse representations of the same single-site coefficients
  /// must produce identical derivative values and validity results.
  #[test]
  fn test_is_zero_branch_optimal_dense_sparse_derivative_parity() {
    let coeffs = vec![0.3, 0.5, 0.1, 0.1];
    let dense = make_dense_contribution(array![[coeffs[0], coeffs[1], coeffs[2], coeffs[3]]]);
    let sparse = make_sparse_contribution(vec![(1.0, coeffs)]);

    assert_eq!(dense.all_sites_valid_at_zero(), sparse.all_sites_valid_at_zero());
    assert_abs_diff_eq!(
      dense.zero_branch_length_derivative(),
      sparse.zero_branch_length_derivative(),
      epsilon = 1e-10
    );
  }

  /// Multi-site parity: 3 sites with different coefficient patterns.
  #[test]
  fn test_is_zero_branch_optimal_dense_sparse_derivative_parity_multi_site() {
    #[rustfmt::skip]
    let dense = make_dense_contribution(array![
      [0.3, 0.5, 0.1, 0.1],
      [0.0, 0.0, 0.8, 0.2],
      [0.1, 0.1, 0.1, 0.7],
    ]);
    let sparse = make_sparse_contribution(vec![
      (1.0, vec![0.3, 0.5, 0.1, 0.1]),
      (1.0, vec![0.0, 0.0, 0.8, 0.2]),
      (1.0, vec![0.1, 0.1, 0.1, 0.7]),
    ]);

    assert_eq!(dense.all_sites_valid_at_zero(), sparse.all_sites_valid_at_zero());
    assert_abs_diff_eq!(
      dense.zero_branch_length_derivative(),
      sparse.zero_branch_length_derivative(),
      epsilon = 1e-10
    );
  }

  /// Create a dense contribution using K80 (Kimura 2-parameter) model.
  ///
  /// K80 has multiple distinct nonzero eigenvalues, so L(t) can be multimodal.
  /// The zero-branch shortcut must not apply for this model.
  fn make_k80_dense_contribution(coefficients: Array2<f64>) -> OptimizationContribution {
    let gtr = k80(K80Params::default()).unwrap();
    OptimizationContribution::Dense(optimize_dense::PartitionContribution::new(coefficients, gtr))
  }

  /// JC69 model is marked unimodal.
  #[test]
  fn test_is_zero_branch_optimal_jc69_is_unimodal() {
    let gtr = jc69(JC69Params::default()).unwrap();
    assert!(gtr.unimodal_branch_likelihood);
  }

  /// F81 model is marked unimodal.
  #[test]
  fn test_is_zero_branch_optimal_f81_is_unimodal() {
    let gtr = f81(F81Params::default()).unwrap();
    assert!(gtr.unimodal_branch_likelihood);
  }

  /// K80 model is NOT marked unimodal (can have two local maxima per
  /// Dinh & Matsen 2017).
  #[test]
  fn test_is_zero_branch_optimal_k80_is_not_unimodal() {
    let gtr = k80(K80Params::default()).unwrap();
    assert!(!gtr.unimodal_branch_likelihood);
  }

  /// All model constructors must classify unimodal_branch_likelihood correctly.
  /// JC69 and F81 have one distinct nonzero eigenvalue (unimodal per Dinh &
  /// Matsen 2017, Corollary 3.1). All others have multiple distinct nonzero
  /// eigenvalues and can exhibit multimodal branch-length likelihood.
  #[rustfmt::skip]
  #[rstest]
  #[case::jc69( jc69(JC69Params::default()).unwrap(),  true)]
  #[case::f81(  f81(F81Params::default()).unwrap(),    true)]
  #[case::k80(  k80(K80Params::default()).unwrap(),    false)]
  #[case::hky85(hky85(HKY85Params::default()).unwrap(), false)]
  #[case::t92(  t92(T92Params::default()).unwrap(),    false)]
  #[case::tn93( tn93(TN93Params::default()).unwrap(),  false)]
  #[case::jtt92(jtt92(Jtt92Params::default()).unwrap(), false)]
  #[trace]
  fn test_is_zero_branch_optimal_model_unimodal_classification(
    #[case] gtr: GTR,
    #[case] expected_unimodal: bool,
  ) {
    assert_eq!(expected_unimodal, gtr.unimodal_branch_likelihood);
  }

  /// K80 contribution with negative derivative at t=0: shortcut returns false
  /// because K80 can be multimodal. The same coefficients under JC69 would
  /// return true.
  #[test]
  fn test_is_zero_branch_optimal_k80_bypasses_shortcut() {
    let k80_contrib = make_k80_dense_contribution(array![[0.0, 1.0, 0.0, 0.0]]);
    // K80 has negative derivative at t=0 (weight on negative eigenvalue)
    assert!(k80_contrib.zero_branch_length_derivative() < 0.0);
    // But shortcut is bypassed because K80 is not proven unimodal
    assert!(!is_zero_branch_optimal(&[k80_contrib]));
  }

  /// Mixed JC69 + K80 contributions: shortcut returns false because at least
  /// one partition uses a non-unimodal model.
  #[test]
  fn test_is_zero_branch_optimal_mixed_unimodal_non_unimodal_returns_false() {
    let jc69_contrib = make_dense_contribution(array![[0.0, 1.0, 0.0, 0.0]]);
    let k80_contrib = make_k80_dense_contribution(array![[0.0, 1.0, 0.0, 0.0]]);
    // JC69 alone would return true
    assert!(is_zero_branch_optimal(&[make_dense_contribution(array![[
      0.0, 1.0, 0.0, 0.0
    ]])]));
    // But mixed with K80, the shortcut is bypassed
    assert!(!is_zero_branch_optimal(&[jc69_contrib, k80_contrib]));
  }

  /// F81 with non-uniform frequencies still uses the shortcut (F81 has one
  /// distinct nonzero eigenvalue regardless of frequency distribution).
  #[test]
  fn test_is_zero_branch_optimal_f81_nonuniform_uses_shortcut() {
    let gtr = f81(F81Params {
      pi: Some(array![0.1, 0.2, 0.3, 0.4]),
      ..F81Params::default()
    })
    .unwrap();
    assert!(gtr.unimodal_branch_likelihood);

    // Coefficients with weight on a negative eigenvalue
    let coefficients = array![[0.0, 1.0, 0.0, 0.0]];
    let contribution = OptimizationContribution::Dense(optimize_dense::PartitionContribution::new(coefficients, gtr));
    assert!(is_zero_branch_optimal(&[contribution]));
  }

  /// K80 sparse contribution also bypasses the shortcut.
  #[test]
  fn test_is_zero_branch_optimal_k80_sparse_bypasses_shortcut() {
    let gtr = k80(K80Params::default()).unwrap();

    let site_contributions = vec![optimize_sparse::SiteContribution {
      multiplicity: 1.0,
      coefficients: array![0.0, 1.0, 0.0, 0.0],
    }];
    let contribution = OptimizationContribution::Sparse(optimize_sparse::PartitionContribution {
      site_contributions,
      gtr,
    });
    // Derivative is negative
    assert!(contribution.zero_branch_length_derivative() < 0.0);
    // But shortcut is bypassed
    assert!(!is_zero_branch_optimal(&[contribution]));
  }

  /// Dinh & Matsen (2017, Section 5) counterexample: K80 with kappa=3 produces
  /// a multimodal branch-length likelihood with three stationary points.
  ///
  /// The construction uses two alignment sites both observed as state A at a
  /// pendant tip, with specific observation probability vectors b^1, b^2 at the
  /// inner node. The per-site likelihood lambda_s(t) = sum_i b^s(i) * P_i(t)
  /// where P_i(t) are K80 transition probabilities from state A.
  ///
  /// K80 kappa=3 eigenvalues: 0, -gamma, -2*gamma where gamma=0.5.
  /// In sorted 4-eigenvalue form: [-1.0, -0.5, -0.5, 0.0].
  ///
  /// The eigenvalue-space coefficients are:
  ///   k_fast  = 0.5 * (b_A - b_T)           [eigenvalue -1.0]
  ///   k_slow  = 0.25 * (b_A + b_T - b_G - b_C)  [eigenvalue -0.5, split across two slots]
  ///   k_const = 0.25                         [eigenvalue 0.0]
  ///
  /// This test verifies that the log-likelihood has three stationary points
  /// (local max, local min, approach to a second maximum at equilibrium),
  /// demonstrating the multimodality that makes the derivative-at-zero
  /// shortcut unreliable for K80.
  #[test]
  fn test_is_zero_branch_optimal_k80_dinh_matsen_multimodal_counterexample() {
    // K80 kappa=3 eigenvalues from the paper (Section 5, eq 5.1)
    let mut gtr = jc69(JC69Params::default()).unwrap();
    gtr.eigvals = array![-1.0, -0.5, -0.5, 0.0];
    gtr.unimodal_branch_likelihood = false;

    // b^1 = [0.24977275, 0.34067358, 0.2051904, 0.20436327] (A, T, G, C)
    // b^2 = [0.25, 0.16087344, 0.29328435, 0.29584221]
    // Eigenvalue-space coefficients derived from the transition probability
    // decomposition (see Section 5, eq 5.2):
    //   k_fast  = 0.5 * (b_A - b_T)
    //   k_slow  = 0.25 * (b_A + b_T - b_G - b_C), split across two -0.5 slots
    //   k_const = 0.25
    #[rustfmt::skip]
    let coefficients = array![
      [-0.04545042, 0.02261158, 0.02261158, 0.25],
      [ 0.04456328, -0.02228164, -0.02228164, 0.25],
    ];
    let contribution = OptimizationContribution::Dense(optimize_dense::PartitionContribution::new(coefficients, gtr));
    let contributions = [contribution];

    // Verify both sites have positive likelihood at t=0
    assert!(contributions[0].all_sites_valid_at_zero());

    // The log-likelihood has three stationary points: a local max near t=0.2,
    // a local min near t=1, and a global max at equilibrium (t -> infinity).
    // This non-monotonic shape -- rise, dip, recovery -- is the hallmark of
    // multimodality that JC69/F81 cannot exhibit.
    let lh = |t: f64| evaluate_mixed_log_lh_only(&contributions, t);

    let log_lh_near_peak = lh(0.2);
    let log_lh_at_dip = lh(1.0);
    let log_lh_at_recovery = lh(5.0);

    // The first peak (near t=0.2) exceeds the dip (near t=1)
    assert!(
      log_lh_near_peak > log_lh_at_dip,
      "K80 counterexample: log_lh near first peak ({log_lh_near_peak}) should exceed dip ({log_lh_at_dip})"
    );
    // Recovery from the dip toward equilibrium
    assert!(
      log_lh_at_recovery > log_lh_at_dip,
      "K80 counterexample: log_lh at recovery ({log_lh_at_recovery}) should exceed dip ({log_lh_at_dip})"
    );

    // The shortcut must return false for K80 regardless of derivative sign
    assert!(!is_zero_branch_optimal(&contributions));
  }
}
