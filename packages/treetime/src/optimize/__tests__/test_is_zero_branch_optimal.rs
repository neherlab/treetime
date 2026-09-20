#[cfg(test)]
mod tests {
  use crate::gtr::get_gtr::{
    F81Params, HKY85Params, JC69Params, Jtt92Params, K80Params, T92Params, TN93Params, f81, hky85, jc69, jtt92, k80,
    t92, tn93,
  };
  use crate::gtr::gtr::GTR;
  use crate::optimize::likelihood::evaluate_mixed_log_lh_only;
  use crate::optimize::zero_boundary::is_zero_branch_optimal;
  use crate::partition::optimize;
  use crate::partition::optimize::contribution::OptimizationContribution;
  use approx::assert_abs_diff_eq;
  use ndarray::{Array2, array};
  use rstest::rstest;

  fn make_dense_contribution(coefficients: Array2<f64>) -> OptimizationContribution {
    let gtr = jc69(JC69Params::default()).unwrap();
    OptimizationContribution::Dense(optimize::dense::PartitionContribution::new(coefficients, gtr))
  }

  fn make_sparse_contribution(sites: Vec<(f64, Vec<f64>)>) -> OptimizationContribution {
    let gtr = jc69(JC69Params::default()).unwrap();

    let site_contributions = sites
      .into_iter()
      .map(|(multiplicity, coeffs)| optimize::sparse::SiteContribution {
        multiplicity,
        coefficients: ndarray::Array1::from_vec(coeffs),
      })
      .collect();

    OptimizationContribution::Sparse(optimize::sparse::PartitionContribution {
      site_contributions,
      gtr,
    })
  }

  #[test]
  fn test_is_zero_branch_optimal_empty_contributions() {
    assert!(!is_zero_branch_optimal(&[]));
  }

  #[test]
  fn test_is_zero_branch_optimal_negative_derivative_returns_true() {
    let contribution = make_dense_contribution(array![[0.0, 1.0, 0.0, 0.0]]);
    assert!(is_zero_branch_optimal(&[contribution]));
  }

  #[test]
  fn test_is_zero_branch_optimal_zero_derivative_returns_false() {
    let contribution = make_dense_contribution(array![[0.0, 0.0, 0.0, 1.0]]);
    assert!(!is_zero_branch_optimal(&[contribution]));
  }

  #[test]
  fn test_is_zero_branch_optimal_scale_invariant_dense() {
    let one_site = make_dense_contribution(array![[0.0, 1.0, 0.0, 0.0]]);
    let result_one = is_zero_branch_optimal(&[one_site]);

    let many_rows: Vec<[f64; 4]> = vec![[0.0, 1.0, 0.0, 0.0]; 100];
    let arr = Array2::from(many_rows);
    let many_sites = make_dense_contribution(arr);
    let result_many = is_zero_branch_optimal(&[many_sites]);

    assert_eq!(
      result_one, result_many,
      "decision must not depend on number of identical sites"
    );
  }

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

  #[test]
  fn test_is_zero_branch_optimal_underflow_resistant() {
    let many_rows: Vec<[f64; 4]> = vec![[0.0625, 0.0625, 0.0625, 0.0625]; 1000];
    let arr = Array2::from(many_rows);
    let contribution = make_dense_contribution(arr);

    assert!(is_zero_branch_optimal(&[contribution]));
  }

  #[test]
  fn test_is_zero_branch_optimal_zero_site_lh_returns_false() {
    let coefficients = array![[0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0],];
    let contribution = make_dense_contribution(coefficients);
    assert!(!is_zero_branch_optimal(&[contribution]));
  }

  #[test]
  fn test_is_zero_branch_optimal_negative_site_lh_returns_false() {
    let coefficients = array![[1.0, -2.0, 0.0, 0.0]];
    let contribution = make_dense_contribution(coefficients);
    assert!(!is_zero_branch_optimal(&[contribution]));
  }

  #[test]
  fn test_is_zero_branch_optimal_nonfinite_site_lh_returns_false() {
    let coefficients = array![[f64::INFINITY, 0.0, 0.0, 0.0]];
    let contribution = make_dense_contribution(coefficients);
    assert!(!is_zero_branch_optimal(&[contribution]));
  }

  #[test]
  fn test_is_zero_branch_optimal_multiple_partitions_negative_derivative() {
    let contrib1 = make_dense_contribution(array![[0.0, 0.5, 0.0, 0.0]]);
    let contrib2 = make_dense_contribution(array![[0.0, 0.5, 0.0, 0.0]]);
    assert!(is_zero_branch_optimal(&[contrib1, contrib2]));
  }

  #[test]
  fn test_is_zero_branch_optimal_small_lh_still_decides() {
    let contrib1 = make_dense_contribution(array![[0.0, 0.05, 0.0, 0.0]]);
    let contrib2 = make_dense_contribution(array![[0.0, 0.05, 0.0, 0.0]]);

    assert!(is_zero_branch_optimal(&[contrib1, contrib2]));
  }

  #[test]
  fn test_is_zero_branch_optimal_mixed_dense_and_sparse() {
    let dense = make_dense_contribution(array![[0.0, 0.5, 0.0, 0.0]]);
    let sparse = make_sparse_contribution(vec![(1.0, vec![0.0, 0.5, 0.0, 0.0])]);
    assert!(is_zero_branch_optimal(&[dense, sparse]));
  }

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

  #[test]
  fn test_is_zero_branch_optimal_sparse_negative_derivative() {
    let contribution = make_sparse_contribution(vec![(1.0, vec![0.0, 1.0, 0.0, 0.0])]);
    assert!(is_zero_branch_optimal(&[contribution]));
  }

  #[test]
  fn test_is_zero_branch_optimal_sparse_zero_derivative() {
    let contribution = make_sparse_contribution(vec![(1.0, vec![0.0, 0.0, 0.0, 1.0])]);
    assert!(!is_zero_branch_optimal(&[contribution]));
  }

  #[test]
  fn test_is_zero_branch_optimal_sparse_multiplicity_preserves_sign() {
    let contribution = make_sparse_contribution(vec![(50.0, vec![0.0, 0.5, 0.0, 0.0])]);
    assert!(is_zero_branch_optimal(&[contribution]));
  }

  #[test]
  fn test_is_zero_branch_optimal_nonfinite_derivative_returns_false() {
    let mut gtr = jc69(JC69Params::default()).unwrap();
    gtr.eigvals = array![-1e300, 0.0, 0.0, 0.0];
    let coefficients = array![[1.0, 0.0, 0.0, -1.0 + f64::EPSILON * 2.0]];
    let contribution = OptimizationContribution::Dense(optimize::dense::PartitionContribution::new(coefficients, gtr));

    assert!(contribution.all_sites_valid_at_zero());
    assert!(!contribution.zero_branch_length_derivative().is_finite());
    assert!(!is_zero_branch_optimal(&[contribution]));
  }

  #[test]
  fn test_is_zero_branch_optimal_derivative_magnitude_dense() {
    let contribution = make_dense_contribution(array![[0.0, 1.0, 0.0, 0.0]]);
    let derivative = contribution.zero_branch_length_derivative();
    assert_abs_diff_eq!(-4.0 / 3.0, derivative, epsilon = 1e-10);
  }

  #[test]
  fn test_is_zero_branch_optimal_derivative_magnitude_sparse() {
    let contribution = make_sparse_contribution(vec![(3.0, vec![0.0, 1.0, 0.0, 0.0])]);
    let derivative = contribution.zero_branch_length_derivative();
    assert_abs_diff_eq!(-4.0, derivative, epsilon = 1e-10);
  }

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

  fn make_k80_dense_contribution(coefficients: Array2<f64>) -> OptimizationContribution {
    let gtr = k80(K80Params::default()).unwrap();
    OptimizationContribution::Dense(optimize::dense::PartitionContribution::new(coefficients, gtr))
  }

  #[test]
  fn test_is_zero_branch_optimal_jc69_is_unimodal() {
    let gtr = jc69(JC69Params::default()).unwrap();
    assert!(gtr.unimodal_branch_likelihood);
  }

  #[test]
  fn test_is_zero_branch_optimal_f81_is_unimodal() {
    let gtr = f81(F81Params::default()).unwrap();
    assert!(gtr.unimodal_branch_likelihood);
  }

  #[test]
  fn test_is_zero_branch_optimal_k80_is_not_unimodal() {
    let gtr = k80(K80Params::default()).unwrap();
    assert!(!gtr.unimodal_branch_likelihood);
  }

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

  #[test]
  fn test_is_zero_branch_optimal_k80_bypasses_shortcut() {
    let k80_contrib = make_k80_dense_contribution(array![[0.0, 1.0, 0.0, 0.0]]);
    assert!(k80_contrib.zero_branch_length_derivative() < 0.0);
    assert!(!is_zero_branch_optimal(&[k80_contrib]));
  }

  #[test]
  fn test_is_zero_branch_optimal_mixed_unimodal_non_unimodal_returns_false() {
    let jc69_contrib = make_dense_contribution(array![[0.0, 1.0, 0.0, 0.0]]);
    let k80_contrib = make_k80_dense_contribution(array![[0.0, 1.0, 0.0, 0.0]]);
    assert!(is_zero_branch_optimal(&[make_dense_contribution(array![[
      0.0, 1.0, 0.0, 0.0
    ]])]));
    assert!(!is_zero_branch_optimal(&[jc69_contrib, k80_contrib]));
  }

  #[test]
  fn test_is_zero_branch_optimal_f81_nonuniform_uses_shortcut() {
    let gtr = f81(F81Params {
      pi: Some(array![0.1, 0.2, 0.3, 0.4]),
      ..F81Params::default()
    })
    .unwrap();
    assert!(gtr.unimodal_branch_likelihood);

    let coefficients = array![[0.0, 1.0, 0.0, 0.0]];
    let contribution = OptimizationContribution::Dense(optimize::dense::PartitionContribution::new(coefficients, gtr));
    assert!(is_zero_branch_optimal(&[contribution]));
  }

  #[test]
  fn test_is_zero_branch_optimal_k80_sparse_bypasses_shortcut() {
    let gtr = k80(K80Params::default()).unwrap();

    let site_contributions = vec![optimize::sparse::SiteContribution {
      multiplicity: 1.0,
      coefficients: array![0.0, 1.0, 0.0, 0.0],
    }];
    let contribution = OptimizationContribution::Sparse(optimize::sparse::PartitionContribution {
      site_contributions,
      gtr,
    });
    assert!(contribution.zero_branch_length_derivative() < 0.0);
    assert!(!is_zero_branch_optimal(&[contribution]));
  }

  #[test]
  fn test_is_zero_branch_optimal_k80_dinh_matsen_multimodal_counterexample() {
    let mut gtr = jc69(JC69Params::default()).unwrap();
    gtr.eigvals = array![-1.0, -0.5, -0.5, 0.0];
    gtr.unimodal_branch_likelihood = false;

    #[rustfmt::skip]
    let coefficients = array![
      [-0.04545042, 0.02261158, 0.02261158, 0.25],
      [ 0.04456328, -0.02228164, -0.02228164, 0.25],
    ];
    let contribution = OptimizationContribution::Dense(optimize::dense::PartitionContribution::new(coefficients, gtr));
    let contributions = [contribution];

    assert!(contributions[0].all_sites_valid_at_zero());

    let lh = |t: f64| {
      evaluate_mixed_log_lh_only(&contributions, t)
        .expect("valid branch length")
        .value()
    };

    let log_lh_near_peak = lh(0.2);
    let log_lh_at_dip = lh(1.0);
    let log_lh_at_recovery = lh(5.0);

    assert!(
      log_lh_near_peak > log_lh_at_dip,
      "K80 counterexample: log_lh near first peak ({log_lh_near_peak}) should exceed dip ({log_lh_at_dip})"
    );
    assert!(
      log_lh_at_recovery > log_lh_at_dip,
      "K80 counterexample: log_lh at recovery ({log_lh_at_recovery}) should exceed dip ({log_lh_at_dip})"
    );

    assert!(!is_zero_branch_optimal(&contributions));
  }
}
