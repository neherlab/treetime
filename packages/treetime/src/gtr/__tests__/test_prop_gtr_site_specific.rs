#[cfg(test)]
mod tests {
  use crate::gtr::gtr::{GTR, GTRParams};
  use crate::gtr::gtr_site_specific::{GTRSiteSpecific, GTRSiteSpecificParams};
  use ndarray::prelude::*;
  use proptest::prelude::*;
  use treetime_utils::array::batched::{matmul_3d, matvec_3d};
  use treetime_utils::{
    pretty_assert_abs_diff_eq, prop_assert_array_abs_diff_eq, prop_assert_array_finite, prop_assert_array_nonneg,
    prop_assert_array_positive, prop_assert_array_upper_bounded,
  };

  mod generators {
    use crate::gtr::__tests__::generators::tests::generators::{arb_pi_nuc, arb_w_nuc};
    use crate::gtr::gtr_site_specific::{GTRSiteSpecific, GTRSiteSpecificParams};
    use ndarray::prelude::*;
    use proptest::prelude::*;

    /// Generate a random site-specific GTR model with the given sequence length.
    pub fn arb_gtr_site_specific(seq_len: usize) -> impl Strategy<Value = GTRSiteSpecific> {
      (
        arb_w_nuc(),
        prop::collection::vec(arb_pi_nuc(), seq_len),
        prop::collection::vec(0.1_f64..5.0, seq_len),
      )
        .prop_map(move |(w, pis, mus)| {
          let mut pi = Array2::zeros((4, seq_len));
          for (a, p) in pis.iter().enumerate() {
            pi.column_mut(a).assign(p);
          }
          let mu = Array1::from_vec(mus);
          GTRSiteSpecific::new(GTRSiteSpecificParams {
            n_states: 4,
            seq_len,
            mu,
            W: Some(w),
            pi,
            approximate: false,
          })
          .expect("GTRSiteSpecific construction should succeed")
        })
    }

    /// Generate a random site-specific GTR model WITH interpolation enabled.
    pub fn arb_gtr_site_specific_approx(seq_len: usize) -> impl Strategy<Value = GTRSiteSpecific> {
      (
        arb_w_nuc(),
        prop::collection::vec(arb_pi_nuc(), seq_len),
        prop::collection::vec(0.1_f64..5.0, seq_len),
      )
        .prop_map(move |(w, pis, mus)| {
          let mut pi = Array2::zeros((4, seq_len));
          for (a, p) in pis.iter().enumerate() {
            pi.column_mut(a).assign(p);
          }
          let mu = Array1::from_vec(mus);
          GTRSiteSpecific::new(GTRSiteSpecificParams {
            n_states: 4,
            seq_len,
            mu,
            W: Some(w),
            pi,
            approximate: true,
          })
          .expect("GTRSiteSpecific construction should succeed")
        })
    }
  }

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    /// Each site's P_a(t) must be column-stochastic: columns sum to 1.
    #[test]
    fn test_prop_gtr_site_specific_expqt_column_stochastic(
      gtr in generators::arb_gtr_site_specific(5),
      t in crate::gtr::__tests__::generators::tests::generators::arb_branch_len(),
    ) {
      let p = gtr.expQt(t).unwrap();
      let (_n, _n2, seq_len) = p.dim();
      let col_sums = p.sum_axis(Axis(0));
      prop_assert_array_abs_diff_eq!(col_sums, Array2::ones((4, seq_len)), epsilon = 1e-10);
    }

    /// P_a(0) = I for all sites.
    #[test]
    fn test_prop_gtr_site_specific_expqt_identity_at_zero(gtr in generators::arb_gtr_site_specific(5)) {
      let p = gtr.expQt(0.0).unwrap();
      let identity_3d = Array3::from_shape_fn((4, 4, 5), |(i, j, _)| if i == j { 1.0 } else { 0.0 });
      prop_assert_array_abs_diff_eq!(p, identity_3d, epsilon = 1e-10);
    }

    /// All entries of P_a(t) must be non-negative (transition probabilities).
    #[test]
    fn test_prop_gtr_site_specific_expqt_nonnegative(
      gtr in generators::arb_gtr_site_specific(5),
      t in crate::gtr::__tests__::generators::tests::generators::arb_branch_len(),
    ) {
      let p = gtr.expQt(t).unwrap();
      prop_assert_array_nonneg!(p, epsilon = 1e-14);
    }

    /// Semigroup: P_a(s+t) = P_a(s) * P_a(t) for all sites.
    #[test]
    fn test_prop_gtr_site_specific_expqt_semigroup(
      gtr in generators::arb_gtr_site_specific(3),
      s in 0.001_f64..1.0,
      t in 0.001_f64..1.0,
    ) {
      let p_s = gtr.expQt(s).unwrap();
      let p_t = gtr.expQt(t).unwrap();
      let p_st = gtr.expQt(s + t).unwrap();

      let product = matmul_3d(&p_s, &p_t);
      prop_assert_array_abs_diff_eq!(product, p_st, epsilon = 1e-8);
    }

    /// As t -> infinity, P_a(t)[i,j] -> pi_a[i] for all j (convergence to equilibrium).
    #[test]
    fn test_prop_gtr_site_specific_expqt_convergence(gtr in generators::arb_gtr_site_specific(3)) {
      let p = gtr.expQt(1000.0).unwrap();
      let expected = Array3::from_shape_fn((4, 4, 3), |(i, _j, a)| gtr.pi[[i, a]]);
      prop_assert_array_abs_diff_eq!(p, expected, epsilon = 1e-6);
    }

    /// propagate_profile produces valid output (non-negative, finite).
    #[test]
    fn test_prop_gtr_site_specific_propagate_profile_valid(
      gtr in generators::arb_gtr_site_specific(5),
      t in crate::gtr::__tests__::generators::tests::generators::arb_branch_len(),
    ) {
      let profile = Array2::from_elem((5, 4), 0.25);
      let result = gtr.propagate_profile(&profile, t, false).unwrap();

      prop_assert_eq!(result.shape(), &[5, 4]);
      prop_assert_array_finite!(result);
      prop_assert_array_nonneg!(result, epsilon = 1e-14);
    }

    /// evolve produces valid output (non-negative, finite).
    #[test]
    fn test_prop_gtr_site_specific_evolve_valid(
      gtr in generators::arb_gtr_site_specific(5),
      t in crate::gtr::__tests__::generators::tests::generators::arb_branch_len(),
    ) {
      let profile = Array2::from_elem((5, 4), 0.25);
      let result = gtr.evolve(&profile, t, false).unwrap();

      prop_assert_eq!(result.shape(), &[5, 4]);
      prop_assert_array_finite!(result);
      prop_assert_array_nonneg!(result, epsilon = 1e-14);
    }

    /// Equilibrium is a fixed point of evolve: evolving pi forward returns pi.
    ///
    /// P(t) is column-stochastic: P @ pi = pi. evolve computes profile @ P^T,
    /// which for profile = pi^T gives pi^T @ P^T = (P @ pi)^T = pi^T.
    #[test]
    fn test_prop_gtr_site_specific_equilibrium_fixed_point(
      gtr in generators::arb_gtr_site_specific(3),
      t in 0.01_f64..10.0,
    ) {
      let profile = gtr.pi.t().to_owned();
      let evolved = gtr.evolve(&profile, t, false).unwrap();
      prop_assert_array_abs_diff_eq!(evolved, profile, epsilon = 1e-8);
    }

    /// Interpolation matches direct computation within 1e-2.
    ///
    /// Linear interpolation on a 61-point grid has inherent accuracy limits
    /// scaling with h^2 * |d^2/dt^2 exp(Qt)|. Observed max error ~1.8e-3
    /// for high-rate sites (mu up to 5.0). Correctness invariants (column
    /// stochastic, non-negative, equilibrium) tested separately at 1e-8 to
    /// 1e-10; expQt_raw validated against v0 at 1e-10 via golden master.
    #[test]
    fn test_prop_gtr_site_specific_interpolation_accuracy(
      gtr in generators::arb_gtr_site_specific_approx(3),
      t in 0.01_f64..2.0,
    ) {
      let p_interp = gtr.expQt(t).unwrap();
      let p_direct = gtr.expQt_raw(t);
      prop_assert_array_abs_diff_eq!(p_interp, p_direct, epsilon = 1e-2);
    }

    /// average_rate is positive for all sites.
    #[test]
    fn test_prop_gtr_site_specific_average_rate_positive(gtr in generators::arb_gtr_site_specific(5)) {
      let rates = gtr.average_rate();
      prop_assert_eq!(rates.len(), 5);
      prop_assert_array_positive!(rates);
    }

    /// All entries of P_a(t) are bounded by 1 (transition probabilities).
    #[test]
    fn test_prop_gtr_site_specific_expqt_bounded(
      gtr in generators::arb_gtr_site_specific(5),
      t in crate::gtr::__tests__::generators::tests::generators::arb_branch_len(),
    ) {
      let p = gtr.expQt(t).unwrap();
      prop_assert_array_upper_bounded!(p, bound = 1.0, epsilon = 1e-14);
    }

    /// Stationary distribution is right eigenvector: P_a(t) @ pi_a = pi_a.
    #[test]
    fn test_prop_gtr_site_specific_stationary_preserved(
      gtr in generators::arb_gtr_site_specific(3),
      t in crate::gtr::__tests__::generators::tests::generators::arb_branch_len(),
    ) {
      let p = gtr.expQt(t).unwrap();
      let pi_evolved = matvec_3d(&p, &gtr.pi);
      prop_assert_array_abs_diff_eq!(pi_evolved, gtr.pi, epsilon = 1e-10);
    }

    /// evolve = profile @ P^T, propagate = profile @ P (per site).
    #[test]
    fn test_prop_gtr_site_specific_evolve_transpose_of_propagate(
      gtr in generators::arb_gtr_site_specific(3),
      t in 0.01_f64..1.0,
    ) {
      let qt = gtr.expQt(t).unwrap();
      let mut profile = Array2::from_elem((3, 4), 0.25);
      // Make profile non-uniform so transpose difference is visible
      profile[[0, 0]] = 0.7;
      profile[[0, 1]] = 0.1;
      profile[[0, 2]] = 0.1;
      profile[[0, 3]] = 0.1;

      let propagated = gtr.propagate_profile(&profile, t, false).unwrap();
      let evolved = gtr.evolve(&profile, t, false).unwrap();

      let expected_prop = helpers::profile_times_transition(&profile, &qt);
      let expected_evol = helpers::profile_times_transition_t(&profile, &qt);
      prop_assert_array_abs_diff_eq!(propagated, expected_prop, epsilon = 1e-10);
      prop_assert_array_abs_diff_eq!(evolved, expected_evol, epsilon = 1e-10);
    }

    /// No NaN or Inf in expQt across random models and branch lengths.
    ///
    /// NaN from 0/0 in normalization, sqrt(negative), or inf-inf in eigendecomposition.
    /// Inf from overflow in exp(large) or division by near-zero in V^{-1}.
    #[test]
    fn test_prop_gtr_site_specific_finite(
      gtr in generators::arb_gtr_site_specific(3),
      t in crate::gtr::__tests__::generators::tests::generators::arb_branch_len(),
    ) {
      let p = gtr.expQt(t).unwrap();
      prop_assert_array_finite!(p);
    }

    /// evolve preserves row sums (probability conservation per site).
    #[test]
    fn test_prop_gtr_site_specific_evolve_preserves_probability(
      gtr in generators::arb_gtr_site_specific(5),
      t in crate::gtr::__tests__::generators::tests::generators::arb_branch_len(),
    ) {
      let profile = Array2::from_elem((5, 4), 0.25);
      let evolved = gtr.evolve(&profile, t, false).unwrap();
      let row_sums = evolved.sum_axis(Axis(1));
      prop_assert_array_abs_diff_eq!(row_sums, Array1::ones(5), epsilon = 1e-10);
    }

    /// Approximate mode: column stochastic.
    #[test]
    fn test_prop_gtr_site_specific_approx_column_stochastic(
      gtr in generators::arb_gtr_site_specific_approx(3),
      t in 0.01_f64..2.0,
    ) {
      let p = gtr.expQt(t).unwrap();
      let col_sums = p.sum_axis(Axis(0));
      prop_assert_array_abs_diff_eq!(col_sums, Array2::ones((4, 3)), epsilon = 1e-8);
    }

    /// Approximate mode: non-negative.
    #[test]
    fn test_prop_gtr_site_specific_approx_nonnegative(
      gtr in generators::arb_gtr_site_specific_approx(3),
      t in 0.01_f64..2.0,
    ) {
      let p = gtr.expQt(t).unwrap();
      prop_assert_array_nonneg!(p, epsilon = 1e-10);
    }

    /// Approximate mode: equilibrium preserved.
    #[test]
    fn test_prop_gtr_site_specific_approx_equilibrium(
      gtr in generators::arb_gtr_site_specific_approx(3),
      t in 0.01_f64..2.0,
    ) {
      let profile = gtr.pi.t().to_owned();
      let evolved = gtr.evolve(&profile, t, false).unwrap();
      prop_assert_array_abs_diff_eq!(evolved, profile, epsilon = 1e-2);
    }
  }

  /// Constructor rejects zero-sum pi column.
  #[test]
  fn test_gtr_site_specific_rejects_zero_pi_column() {
    let pi = array![[0.25, 0.0], [0.25, 0.0], [0.25, 0.0], [0.25, 0.0]];
    let result = GTRSiteSpecific::new(GTRSiteSpecificParams {
      n_states: 4,
      seq_len: 2,
      mu: array![1.0, 1.0],
      W: None,
      pi,
      approximate: false,
    });
    assert!(result.is_err(), "Should reject zero-sum pi column");
  }

  /// Constructor rejects negative mu.
  #[test]
  fn test_gtr_site_specific_rejects_negative_mu() {
    let pi = array![[0.25, 0.25], [0.25, 0.25], [0.25, 0.25], [0.25, 0.25]];
    let result = GTRSiteSpecific::new(GTRSiteSpecificParams {
      n_states: 4,
      seq_len: 2,
      mu: array![1.0, -0.5],
      W: None,
      pi,
      approximate: false,
    });
    assert!(result.is_err(), "Should reject negative mu");
  }

  /// Constructor rejects dimension mismatch.
  #[test]
  fn test_gtr_site_specific_rejects_dimension_mismatch() {
    let pi = array![[0.25, 0.25], [0.25, 0.25], [0.25, 0.25], [0.25, 0.25]];
    let result = GTRSiteSpecific::new(GTRSiteSpecificParams {
      n_states: 4,
      seq_len: 2,
      mu: array![1.0, 1.0, 1.0], // 3 != seq_len=2
      W: None,
      pi,
      approximate: false,
    });
    assert!(result.is_err(), "Should reject mu length mismatch");
  }

  /// Verify that different sites produce different transition matrices.
  #[test]
  fn test_gtr_site_specific_different_sites_differ() {
    let gtr = helpers::two_site_model();
    let p = gtr.expQt(1.0).unwrap();
    let p0 = p.slice(s![.., .., 0]);
    let p1 = p.slice(s![.., .., 1]);

    let max_diff = (&p0 - &p1).mapv(f64::abs).iter().copied().fold(0.0_f64, f64::max);
    assert!(
      max_diff > 0.01,
      "Sites with different pi should produce different P(t), max diff = {max_diff}"
    );
  }

  /// When all sites have the same pi, site-specific model matches standard GTR.
  #[test]
  fn test_gtr_site_specific_uniform_matches_standard() {
    let pi = array![0.1, 0.2, 0.3, 0.4];
    let W = {
      let mut w = Array2::<f64>::ones([4, 4]);
      w[[0, 2]] = 3.0;
      w[[2, 0]] = 3.0;
      w[[1, 3]] = 3.0;
      w[[3, 1]] = 3.0;
      w
    };

    let standard = GTR::new(GTRParams {
      n_states: 4,
      mu: 1.0,
      W: Some(W.clone()),
      pi: pi.clone(),
    })
    .unwrap();

    let seq_len = 3;
    let mut pi_2d = Array2::zeros((4, seq_len));
    for a in 0..seq_len {
      pi_2d.column_mut(a).assign(&pi);
    }

    let site_specific = GTRSiteSpecific::new(GTRSiteSpecificParams {
      n_states: 4,
      seq_len,
      mu: Array1::ones(seq_len),
      W: Some(W),
      pi: pi_2d,
      approximate: false,
    })
    .unwrap();

    let t = 0.5;
    let p_standard = standard.expQt(t);
    let p_ss = site_specific.expQt(t).unwrap();

    let expected_3d = Array3::from_shape_fn((4, 4, seq_len), |(i, j, _)| p_standard[[i, j]]);
    pretty_assert_abs_diff_eq!(p_ss, expected_3d, epsilon = 1e-10);
  }

  mod helpers {
    use super::*;

    /// Per-site `profile @ P`: `result[a,j] = sum_i profile[a,i] * P[i,j,a]`.
    ///
    /// Profile shape `(sites, states)`, transition tensor shape `(states, states, sites)`.
    pub fn profile_times_transition(profile: &Array2<f64>, p: &Array3<f64>) -> Array2<f64> {
      let n_sites = profile.nrows();
      let mut result = Array2::zeros(profile.raw_dim());
      for a in 0..n_sites {
        result.row_mut(a).assign(&profile.row(a).dot(&p.slice(s![.., .., a])));
      }
      result
    }

    /// Per-site `profile @ P^T`: `result[a,i] = sum_j profile[a,j] * P[i,j,a]`.
    ///
    /// Profile shape `(sites, states)`, transition tensor shape `(states, states, sites)`.
    pub fn profile_times_transition_t(profile: &Array2<f64>, p: &Array3<f64>) -> Array2<f64> {
      let n_sites = profile.nrows();
      let mut result = Array2::zeros(profile.raw_dim());
      for a in 0..n_sites {
        result
          .row_mut(a)
          .assign(&profile.row(a).dot(&p.slice(s![.., .., a]).t()));
      }
      result
    }

    /// Construct a simple 2-site model where site 0 has uniform pi and site 1 has skewed pi.
    pub fn two_site_model() -> GTRSiteSpecific {
      let pi = array![[0.25, 0.1], [0.25, 0.2], [0.25, 0.3], [0.25, 0.4]];
      GTRSiteSpecific::new(GTRSiteSpecificParams {
        n_states: 4,
        seq_len: 2,
        mu: array![1.0, 1.0],
        W: None,
        pi,
        approximate: false,
      })
      .expect("Construction should succeed")
    }
  }
}
