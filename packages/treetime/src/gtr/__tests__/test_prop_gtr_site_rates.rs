#[cfg(test)]
mod tests {
  use crate::gtr::__tests__::generators::tests::generators::{arb_branch_len, arb_gtr_nuc, arb_profile_nuc};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use approx::assert_abs_diff_eq;
  use ndarray::{Array1, Array2};
  use proptest::prelude::*;
  use treetime_utils::pretty_assert_abs_diff_eq;

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    /// expQt_with_rate(t, 1.0) == expQt(t)
    #[test]
    fn test_prop_expqt_with_rate_one_equals_expqt(gtr in arb_gtr_nuc(), t in arb_branch_len()) {
      let p_uniform = gtr.expQt(t);
      let p_rate_one = gtr.expQt_with_rate(t, 1.0);
      pretty_assert_abs_diff_eq!(p_uniform, p_rate_one, epsilon = 1e-14);
    }

    /// expQt_with_rate(t, r) is column-stochastic (columns sum to 1)
    #[test]
    fn test_prop_expqt_with_rate_column_stochastic(
      gtr in arb_gtr_nuc(),
      t in arb_branch_len(),
      rate in 0.1_f64..5.0,
    ) {
      let p = gtr.expQt_with_rate(t, rate);
      for j in 0..p.ncols() {
        let col_sum = p.column(j).sum();
        prop_assert!(
          (col_sum - 1.0).abs() < 1e-10,
          "column {j} sum = {col_sum}, expected 1.0"
        );
      }
    }

    /// expQt_with_rate(t, r) has non-negative entries
    #[test]
    fn test_prop_expqt_with_rate_non_negative(
      gtr in arb_gtr_nuc(),
      t in arb_branch_len(),
      rate in 0.1_f64..5.0,
    ) {
      let p = gtr.expQt_with_rate(t, rate);
      for &val in &p {
        prop_assert!(val >= 0.0, "P(t) entry must be >= 0, got {val}");
      }
    }

    /// Semigroup: P(r*t1) * P(r*t2) = P(r*(t1+t2)) when rates are uniform
    #[test]
    fn test_prop_expqt_with_rate_semigroup(
      gtr in arb_gtr_nuc(),
      s in 0.001_f64..1.0,
      t in 0.001_f64..1.0,
      rate in 0.5_f64..3.0,
    ) {
      let p_s = gtr.expQt_with_rate(s, rate);
      let p_t = gtr.expQt_with_rate(t, rate);
      let p_st = gtr.expQt_with_rate(s + t, rate);
      let product = p_s.dot(&p_t);
      pretty_assert_abs_diff_eq!(product, p_st, epsilon = 1e-10);
    }

    /// propagate_profile with uniform site_rates == propagate_profile without site_rates
    #[test]
    fn test_prop_propagate_uniform_rates_equals_scalar(
      gtr in arb_gtr_nuc(),
      profile in arb_profile_nuc(10),
      t in arb_branch_len(),
    ) {
      let result_scalar = gtr.propagate_profile(&profile, t, false);

      let mut gtr_per_site = gtr;
      gtr_per_site.set_site_rates(Array1::ones(10));
      let result_per_site = gtr_per_site.propagate_profile(&profile, t, false);

      pretty_assert_abs_diff_eq!(result_scalar, result_per_site, epsilon = 1e-12);
    }

    /// evolve with uniform site_rates == evolve without site_rates
    #[test]
    fn test_prop_evolve_uniform_rates_equals_scalar(
      gtr in arb_gtr_nuc(),
      profile in arb_profile_nuc(10),
      t in arb_branch_len(),
    ) {
      let result_scalar = gtr.evolve(&profile, t, false);

      let mut gtr_per_site = gtr;
      gtr_per_site.set_site_rates(Array1::ones(10));
      let result_per_site = gtr_per_site.evolve(&profile, t, false);

      pretty_assert_abs_diff_eq!(result_scalar, result_per_site, epsilon = 1e-12);
    }

    /// Per-site propagate_profile produces row-wise results matching individual expQt_with_rate
    #[test]
    fn test_prop_propagate_per_site_matches_individual(
      gtr in arb_gtr_nuc(),
      profile in arb_profile_nuc(5),
      t in 0.001_f64..1.0,
      rates in prop::collection::vec(0.2_f64..4.0, 5),
    ) {
      let site_rates = Array1::from_vec(rates.clone());
      let mut gtr_per_site = gtr.clone();
      gtr_per_site.set_site_rates(site_rates);

      let result = gtr_per_site.propagate_profile(&profile, t, false);

      // Compare each row against individual expQt_with_rate computation
      for (l, rate) in rates.iter().enumerate() {
        let exp_qt = gtr.expQt_with_rate(t, *rate);
        let row = profile.row(l);
        let expected = row.dot(&exp_qt);
        for k in 0..4 {
          prop_assert!(
            (result[[l, k]] - expected[k]).abs() < 1e-12,
            "row {l}, col {k}: result={} expected={}", result[[l, k]], expected[k]
          );
        }
      }
    }

    /// Per-site evolve produces row-wise results matching individual expQt_with_rate
    #[test]
    fn test_prop_evolve_per_site_matches_individual(
      gtr in arb_gtr_nuc(),
      profile in arb_profile_nuc(5),
      t in 0.001_f64..1.0,
      rates in prop::collection::vec(0.2_f64..4.0, 5),
    ) {
      let site_rates = Array1::from_vec(rates.clone());
      let mut gtr_per_site = gtr.clone();
      gtr_per_site.set_site_rates(site_rates);

      let result = gtr_per_site.evolve(&profile, t, false);

      for (l, rate) in rates.iter().enumerate() {
        let exp_qt = gtr.expQt_with_rate(t, *rate);
        let row = profile.row(l);
        let expected = row.dot(&exp_qt.t());
        for k in 0..4 {
          prop_assert!(
            (result[[l, k]] - expected[k]).abs() < 1e-12,
            "row {l}, col {k}: result={} expected={}", result[[l, k]], expected[k]
          );
        }
      }
    }

    /// Per-site propagate at t=0 is identity (rates don't matter)
    #[test]
    fn test_prop_propagate_per_site_identity_at_zero(
      gtr in arb_gtr_nuc(),
      profile in arb_profile_nuc(5),
      rates in prop::collection::vec(0.1_f64..10.0, 5),
    ) {
      let mut gtr_per_site = gtr;
      gtr_per_site.set_site_rates(Array1::from_vec(rates));

      let result = gtr_per_site.propagate_profile(&profile, 0.0, false);
      pretty_assert_abs_diff_eq!(result, profile, epsilon = 1e-14);
    }

    /// Large t: evolve with per-site rates converges to equilibrium (each row -> pi)
    #[test]
    fn test_prop_evolve_per_site_equilibrium(
      gtr in arb_gtr_nuc(),
      rates in prop::collection::vec(0.5_f64..3.0, 5),
    ) {
      let min_rate = rates.iter().copied().fold(f64::INFINITY, f64::min);
      let mut gtr_per_site = gtr.clone();
      gtr_per_site.set_site_rates(Array1::from_vec(rates));
      // Scale t by 1/min_rate to ensure even the slowest site converges
      let t = 1000.0 / (gtr.mu.max(0.001) * min_rate);

      let profile = Array2::from_shape_fn((5, 4), |(_l, k)| if k == 0 { 1.0 } else { 0.0 });
      let result = gtr_per_site.evolve(&profile, t, false);

      for l in 0..5 {
        let row = result.row(l);
        for k in 0..4 {
          prop_assert!(
            (row[k] - gtr.pi[k]).abs() < 1e-4,
            "row {l}, col {k}: result={} expected pi={}", row[k], gtr.pi[k]
          );
        }
      }
    }
  }

  /// Higher site rate produces more divergence from initial state
  #[test]
  fn test_higher_rate_more_divergence() {
    let gtr = jc69(JC69Params::default()).unwrap();
    let t = 0.5;

    let p_slow = gtr.expQt_with_rate(t, 0.1);
    let p_fast = gtr.expQt_with_rate(t, 5.0);

    // Diagonal element (staying in same state) should be higher for slow rate
    assert!(
      p_slow[[0, 0]] > p_fast[[0, 0]],
      "slow rate should preserve state more: slow={} fast={}",
      p_slow[[0, 0]],
      p_fast[[0, 0]]
    );
  }

  /// Rate 0 means no evolution (identity matrix)
  #[test]
  fn test_rate_zero_is_identity() {
    let gtr = jc69(JC69Params::default()).unwrap();
    let p = gtr.expQt_with_rate(1.0, 0.0);
    pretty_assert_abs_diff_eq!(p, Array2::eye(4), epsilon = 1e-14);
  }

  /// set_site_rates / clear_site_rates / has_site_rates lifecycle
  #[test]
  fn test_site_rates_lifecycle() {
    let mut gtr = jc69(JC69Params::default()).unwrap();
    assert!(!gtr.has_site_rates());

    gtr.set_site_rates(Array1::from_vec(vec![0.5, 1.0, 1.5, 2.0]));
    assert!(gtr.has_site_rates());
    assert_abs_diff_eq!(gtr.site_rates.as_ref().unwrap()[2], 1.5, epsilon = 1e-14);

    gtr.clear_site_rates();
    assert!(!gtr.has_site_rates());
  }
}
