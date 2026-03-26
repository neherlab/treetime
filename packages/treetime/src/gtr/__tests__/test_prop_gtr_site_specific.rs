#[cfg(test)]
mod tests {
  use crate::gtr::__tests__::generators::tests::generators::{arb_branch_len, arb_pi_nuc, arb_w_nuc};
  use crate::gtr::gtr_site_specific::{GTRSiteSpecific, GTRSiteSpecificParams};
  use ndarray::prelude::*;
  use proptest::prelude::*;

  /// Generate a random site-specific GTR model with the given sequence length.
  fn arb_gtr_site_specific(seq_len: usize) -> impl Strategy<Value = GTRSiteSpecific> {
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
  fn arb_gtr_site_specific_approx(seq_len: usize) -> impl Strategy<Value = GTRSiteSpecific> {
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

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    /// Each site's P_a(t) must be column-stochastic: columns sum to 1.
    #[test]
    fn test_prop_gtr_site_specific_expqt_column_stochastic(
      gtr in arb_gtr_site_specific(5),
      t in arb_branch_len(),
    ) {
      let p = gtr.expQt(t);
      let (_n, _n2, seq_len) = p.dim();
      for a in 0..seq_len {
        let p_a = p.slice(s![.., .., a]);
        for j in 0..4 {
          let col_sum: f64 = p_a.column(j).sum();
          prop_assert!(
            (col_sum - 1.0).abs() < 1e-10,
            "P_a(t={t})[site={a}] column {j} sum = {col_sum}, expected 1.0"
          );
        }
      }
    }

    /// P_a(0) = I for all sites.
    #[test]
    fn test_prop_gtr_site_specific_expqt_identity_at_zero(gtr in arb_gtr_site_specific(5)) {
      let p = gtr.expQt(0.0);
      let identity = Array2::<f64>::eye(4);
      for a in 0..5 {
        let p_a = p.slice(s![.., .., a]);
        for i in 0..4 {
          for j in 0..4 {
            prop_assert!(
              (p_a[[i, j]] - identity[[i, j]]).abs() < 1e-10,
              "P(0)[{i},{j},site={a}] = {}, expected {}",
              p_a[[i, j]],
              identity[[i, j]]
            );
          }
        }
      }
    }

    /// All entries of P_a(t) must be non-negative (transition probabilities).
    #[test]
    fn test_prop_gtr_site_specific_expqt_nonnegative(
      gtr in arb_gtr_site_specific(5),
      t in arb_branch_len(),
    ) {
      let p = gtr.expQt(t);
      for a in 0..5 {
        for i in 0..4 {
          for j in 0..4 {
            prop_assert!(
              p[[i, j, a]] >= -1e-14,
              "P(t={t})[{i},{j},site={a}] = {} is negative",
              p[[i, j, a]]
            );
          }
        }
      }
    }

    /// Semigroup: P_a(s+t) = P_a(s) * P_a(t) for all sites.
    #[test]
    fn test_prop_gtr_site_specific_expqt_semigroup(
      gtr in arb_gtr_site_specific(3),
      s in 0.001_f64..1.0,
      t in 0.001_f64..1.0,
    ) {
      let p_s = gtr.expQt(s);
      let p_t = gtr.expQt(t);
      let p_st = gtr.expQt(s + t);

      for a in 0..3 {
        let ps_a: Array2<f64> = p_s.slice(s![.., .., a]).to_owned();
        let pt_a: Array2<f64> = p_t.slice(s![.., .., a]).to_owned();
        let pst_a = p_st.slice(s![.., .., a]);
        let product = ps_a.dot(&pt_a);
        for i in 0..4 {
          for j in 0..4 {
            prop_assert!(
              (product[[i, j]] - pst_a[[i, j]]).abs() < 1e-8,
              "Semigroup violation at [{i},{j},site={a}]: P(s)*P(t) = {}, P(s+t) = {}",
              product[[i, j]],
              pst_a[[i, j]]
            );
          }
        }
      }
    }

    /// As t -> infinity, P_a(t)[i,j] -> pi_a[i] for all j (convergence to equilibrium).
    #[test]
    fn test_prop_gtr_site_specific_expqt_convergence(gtr in arb_gtr_site_specific(3)) {
      let p = gtr.expQt(1000.0);
      for a in 0..3 {
        let p_a = p.slice(s![.., .., a]);
        let pi_a = gtr.pi.column(a);
        for j in 0..4 {
          for i in 0..4 {
            prop_assert!(
              (p_a[[i, j]] - pi_a[i]).abs() < 1e-6,
              "Convergence: P(1000)[{i},{j},site={a}] = {}, pi[{i}] = {}",
              p_a[[i, j]],
              pi_a[i]
            );
          }
        }
      }
    }

    /// propagate_profile produces valid output (non-negative, finite).
    #[test]
    fn test_prop_gtr_site_specific_propagate_profile_valid(
      gtr in arb_gtr_site_specific(5),
      t in arb_branch_len(),
    ) {
      // Create a simple profile: each site has a uniform distribution
      let profile = Array2::from_elem((5, 4), 0.25);
      let result = gtr.propagate_profile(&profile, t, false);

      prop_assert_eq!(result.shape(), &[5, 4]);
      for a in 0..5 {
        for i in 0..4 {
          let v = result[[a, i]];
          prop_assert!(v.is_finite(), "propagate_profile[{a},{i}] is not finite: {v}");
          prop_assert!(v >= -1e-14, "propagate_profile[{a},{i}] is negative: {v}");
        }
      }
    }

    /// evolve produces valid output (non-negative, finite).
    #[test]
    fn test_prop_gtr_site_specific_evolve_valid(
      gtr in arb_gtr_site_specific(5),
      t in arb_branch_len(),
    ) {
      let profile = Array2::from_elem((5, 4), 0.25);
      let result = gtr.evolve(&profile, t, false);

      prop_assert_eq!(result.shape(), &[5, 4]);
      for a in 0..5 {
        for i in 0..4 {
          let v = result[[a, i]];
          prop_assert!(v.is_finite(), "evolve[{a},{i}] is not finite: {v}");
          prop_assert!(v >= -1e-14, "evolve[{a},{i}] is negative: {v}");
        }
      }
    }

    /// Equilibrium is a fixed point of evolve: evolving pi forward returns pi.
    ///
    /// P(t) is column-stochastic: P @ pi = pi. evolve computes profile @ P^T,
    /// which for profile = pi^T gives pi^T @ P^T = (P @ pi)^T = pi^T.
    #[test]
    fn test_prop_gtr_site_specific_equilibrium_fixed_point(
      gtr in arb_gtr_site_specific(3),
      t in 0.01_f64..10.0,
    ) {
      // Build profile where each site's row is that site's pi
      let mut profile = Array2::zeros((3, 4));
      for a in 0..3 {
        profile.row_mut(a).assign(&gtr.pi.column(a));
      }

      let evolved = gtr.evolve(&profile, t, false);
      for a in 0..3 {
        for i in 0..4 {
          prop_assert!(
            (evolved[[a, i]] - profile[[a, i]]).abs() < 1e-8,
            "Equilibrium not fixed at [{a},{i}]: evolved = {}, expected = {}",
            evolved[[a, i]],
            profile[[a, i]]
          );
        }
      }
    }

    /// Interpolation matches direct computation within tolerance.
    ///
    /// Linear interpolation of the matrix exponential introduces error proportional
    /// to the grid spacing squared and the second derivative of exp(Qt). For models
    /// with high per-site rates, the function changes rapidly and 5e-3 tolerance
    /// is appropriate for the 61-point non-uniform grid used by the interpolator.
    #[test]
    fn test_prop_gtr_site_specific_interpolation_accuracy(
      gtr in arb_gtr_site_specific_approx(3),
      t in 0.01_f64..2.0,
    ) {
      let p_interp = gtr.expQt(t);
      let p_direct = gtr.expQt_raw(t);

      for a in 0..3 {
        for i in 0..4 {
          for j in 0..4 {
            prop_assert!(
              (p_interp[[i, j, a]] - p_direct[[i, j, a]]).abs() < 5e-3,
              "Interpolation error at [{i},{j},site={a}], t={t}: interp = {}, direct = {}",
              p_interp[[i, j, a]],
              p_direct[[i, j, a]]
            );
          }
        }
      }
    }

    /// average_rate is positive for all sites.
    #[test]
    fn test_prop_gtr_site_specific_average_rate_positive(gtr in arb_gtr_site_specific(5)) {
      let rates = gtr.average_rate();
      prop_assert_eq!(rates.len(), 5);
      for a in 0..5 {
        prop_assert!(rates[a] > 0.0, "average_rate[{a}] = {} is not positive", rates[a]);
      }
    }
  }

  mod helpers {
    use super::*;

    /// Construct a simple 2-site model where site 0 has uniform pi and site 1 has skewed pi.
    /// This enables testing that per-site eigendecomposition produces genuinely different
    /// P(t) matrices for different sites.
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

  /// Verify that different sites produce different transition matrices.
  #[test]
  fn test_gtr_site_specific_different_sites_differ() {
    let gtr = helpers::two_site_model();
    let p = gtr.expQt(1.0);
    let p0 = p.slice(s![.., .., 0]);
    let p1 = p.slice(s![.., .., 1]);

    // Sites should differ because they have different pi
    let max_diff = (&p0 - &p1).mapv(f64::abs).iter().copied().fold(0.0_f64, f64::max);
    assert!(
      max_diff > 0.01,
      "Sites with different pi should produce different P(t), max diff = {max_diff}"
    );
  }

  /// When all sites have the same pi, site-specific model matches standard GTR.
  #[test]
  fn test_gtr_site_specific_uniform_matches_standard() {
    use crate::gtr::gtr::{GTR, GTRParams};

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
    let p_ss = site_specific.expQt(t);

    for a in 0..seq_len {
      let p_a = p_ss.slice(s![.., .., a]);
      for i in 0..4 {
        for j in 0..4 {
          assert!(
            (p_a[[i, j]] - p_standard[[i, j]]).abs() < 1e-10,
            "Site {a}: P_ss[{i},{j}] = {}, P_std[{i},{j}] = {}, diff = {}",
            p_a[[i, j]],
            p_standard[[i, j]],
            (p_a[[i, j]] - p_standard[[i, j]]).abs()
          );
        }
      }
    }
  }
}
