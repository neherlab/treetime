#[cfg(test)]
mod tests {
  use crate::gtr::gtr::{GTR, GTRParams};
  use crate::gtr::gtr_site_specific::{GTRSiteSpecific, GTRSiteSpecificParams};
  use ndarray::prelude::*;
  use proptest::prelude::*;

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
    fn test_prop_gtr_site_specific_expqt_identity_at_zero(gtr in generators::arb_gtr_site_specific(5)) {
      let p = gtr.expQt(0.0).unwrap();
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
      gtr in generators::arb_gtr_site_specific(5),
      t in crate::gtr::__tests__::generators::tests::generators::arb_branch_len(),
    ) {
      let p = gtr.expQt(t).unwrap();
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
      gtr in generators::arb_gtr_site_specific(3),
      s in 0.001_f64..1.0,
      t in 0.001_f64..1.0,
    ) {
      let p_s = gtr.expQt(s).unwrap();
      let p_t = gtr.expQt(t).unwrap();
      let p_st = gtr.expQt(s + t).unwrap();

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
    fn test_prop_gtr_site_specific_expqt_convergence(gtr in generators::arb_gtr_site_specific(3)) {
      let p = gtr.expQt(1000.0).unwrap();
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
      gtr in generators::arb_gtr_site_specific(5),
      t in crate::gtr::__tests__::generators::tests::generators::arb_branch_len(),
    ) {
      let profile = Array2::from_elem((5, 4), 0.25);
      let result = gtr.propagate_profile(&profile, t, false).unwrap();

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
      gtr in generators::arb_gtr_site_specific(5),
      t in crate::gtr::__tests__::generators::tests::generators::arb_branch_len(),
    ) {
      let profile = Array2::from_elem((5, 4), 0.25);
      let result = gtr.evolve(&profile, t, false).unwrap();

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
      gtr in generators::arb_gtr_site_specific(3),
      t in 0.01_f64..10.0,
    ) {
      let mut profile = Array2::zeros((3, 4));
      for a in 0..3 {
        profile.row_mut(a).assign(&gtr.pi.column(a));
      }

      let evolved = gtr.evolve(&profile, t, false).unwrap();
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

      for a in 0..3 {
        for i in 0..4 {
          for j in 0..4 {
            prop_assert!(
              (p_interp[[i, j, a]] - p_direct[[i, j, a]]).abs() < 1e-2,
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
    fn test_prop_gtr_site_specific_average_rate_positive(gtr in generators::arb_gtr_site_specific(5)) {
      let rates = gtr.average_rate();
      prop_assert_eq!(rates.len(), 5);
      for a in 0..5 {
        prop_assert!(rates[a] > 0.0, "average_rate[{a}] = {} is not positive", rates[a]);
      }
    }

    /// All entries of P_a(t) are bounded by 1 (transition probabilities).
    #[test]
    fn test_prop_gtr_site_specific_expqt_bounded(
      gtr in generators::arb_gtr_site_specific(5),
      t in crate::gtr::__tests__::generators::tests::generators::arb_branch_len(),
    ) {
      let p = gtr.expQt(t).unwrap();
      for a in 0..5 {
        for i in 0..4 {
          for j in 0..4 {
            prop_assert!(
              p[[i, j, a]] <= 1.0 + 1e-14,
              "P(t={t})[{i},{j},site={a}] = {} exceeds 1.0",
              p[[i, j, a]]
            );
          }
        }
      }
    }

    /// Stationary distribution is right eigenvector: P_a(t) @ pi_a = pi_a.
    #[test]
    fn test_prop_gtr_site_specific_stationary_preserved(
      gtr in generators::arb_gtr_site_specific(3),
      t in crate::gtr::__tests__::generators::tests::generators::arb_branch_len(),
    ) {
      let p = gtr.expQt(t).unwrap();
      for a in 0..3 {
        let p_a: Array2<f64> = p.slice(s![.., .., a]).to_owned();
        let pi_a = gtr.pi.column(a).to_owned();
        let pi_evolved = p_a.dot(&pi_a);
        for i in 0..4 {
          prop_assert!(
            (pi_evolved[i] - pi_a[i]).abs() < 1e-10,
            "Stationary not preserved at [{i},site={a}]: P@pi = {}, pi = {}",
            pi_evolved[i],
            pi_a[i]
          );
        }
      }
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

      for a in 0..3 {
        let qt_a = qt.slice(s![.., .., a]);
        let expected_prop = profile.row(a).dot(&qt_a);
        let expected_evol = profile.row(a).dot(&qt_a.t());
        for i in 0..4 {
          prop_assert!(
            (propagated[[a, i]] - expected_prop[i]).abs() < 1e-10,
            "propagate mismatch at [{a},{i}]"
          );
          prop_assert!(
            (evolved[[a, i]] - expected_evol[i]).abs() < 1e-10,
            "evolve mismatch at [{a},{i}]"
          );
        }
      }
    }

    /// No NaN in expQt across random models and branch lengths.
    #[test]
    fn test_prop_gtr_site_specific_no_nan(
      gtr in generators::arb_gtr_site_specific(3),
      t in crate::gtr::__tests__::generators::tests::generators::arb_branch_len(),
    ) {
      let p = gtr.expQt(t).unwrap();
      for ((i, j, a), &v) in p.indexed_iter() {
        prop_assert!(!v.is_nan(), "P(t={t})[{i},{j},site={a}] is NaN");
      }
    }

    /// No Inf in expQt across random models and branch lengths.
    #[test]
    fn test_prop_gtr_site_specific_no_inf(
      gtr in generators::arb_gtr_site_specific(3),
      t in crate::gtr::__tests__::generators::tests::generators::arb_branch_len(),
    ) {
      let p = gtr.expQt(t).unwrap();
      for ((i, j, a), &v) in p.indexed_iter() {
        prop_assert!(!v.is_infinite(), "P(t={t})[{i},{j},site={a}] is Inf");
      }
    }

    /// evolve preserves row sums (probability conservation per site).
    #[test]
    fn test_prop_gtr_site_specific_evolve_preserves_probability(
      gtr in generators::arb_gtr_site_specific(5),
      t in crate::gtr::__tests__::generators::tests::generators::arb_branch_len(),
    ) {
      let profile = Array2::from_elem((5, 4), 0.25);
      let evolved = gtr.evolve(&profile, t, false).unwrap();
      for a in 0..5 {
        let row_sum: f64 = evolved.row(a).sum();
        prop_assert!(
          (row_sum - 1.0).abs() < 1e-10,
          "evolve row {a} sum = {row_sum}, expected 1.0"
        );
      }
    }

    /// Approximate mode: column stochastic.
    #[test]
    fn test_prop_gtr_site_specific_approx_column_stochastic(
      gtr in generators::arb_gtr_site_specific_approx(3),
      t in 0.01_f64..2.0,
    ) {
      let p = gtr.expQt(t).unwrap();
      for a in 0..3 {
        let p_a = p.slice(s![.., .., a]);
        for j in 0..4 {
          let col_sum: f64 = p_a.column(j).sum();
          prop_assert!(
            (col_sum - 1.0).abs() < 1e-8,
            "Approx P_a(t={t})[site={a}] column {j} sum = {col_sum}, expected 1.0"
          );
        }
      }
    }

    /// Approximate mode: non-negative.
    #[test]
    fn test_prop_gtr_site_specific_approx_nonnegative(
      gtr in generators::arb_gtr_site_specific_approx(3),
      t in 0.01_f64..2.0,
    ) {
      let p = gtr.expQt(t).unwrap();
      for ((i, j, a), &v) in p.indexed_iter() {
        prop_assert!(v >= -1e-10, "Approx P(t={t})[{i},{j},site={a}] = {v} is negative");
      }
    }

    /// Approximate mode: equilibrium preserved.
    #[test]
    fn test_prop_gtr_site_specific_approx_equilibrium(
      gtr in generators::arb_gtr_site_specific_approx(3),
      t in 0.01_f64..2.0,
    ) {
      let mut profile = Array2::zeros((3, 4));
      for a in 0..3 {
        profile.row_mut(a).assign(&gtr.pi.column(a));
      }
      let evolved = gtr.evolve(&profile, t, false).unwrap();
      for a in 0..3 {
        for i in 0..4 {
          prop_assert!(
            (evolved[[a, i]] - profile[[a, i]]).abs() < 1e-2,
            "Approx equilibrium not preserved at [{a},{i}]: {}, expected {}",
            evolved[[a, i]],
            profile[[a, i]]
          );
        }
      }
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

  mod helpers {
    use super::*;

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
