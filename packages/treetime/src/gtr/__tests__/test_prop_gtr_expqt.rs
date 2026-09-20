#[cfg(test)]
mod tests {
  use crate::gtr::__tests__::generators::tests::generators::{arb_branch_len, arb_gtr_nuc, arb_profile_nuc};
  use crate::gtr::__tests__::prop_support::{prop_assert_columns_sum_to, prop_assert_rows_sum_to};
  use ndarray::Array2;
  use proptest::prelude::*;
  use treetime_utils::{prop_assert_array_abs_diff_eq, prop_assert_array_nonneg, prop_assert_array_upper_bounded};

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    #[test]
    fn test_prop_gtr_expqt_stochastic_columns(gtr in arb_gtr_nuc(), t in arb_branch_len()) {
      let p = gtr.expQt(t);
      prop_assert_columns_sum_to(&p, 1.0, 1e-10)?;
    }

    #[test]
    fn test_prop_gtr_expqt_nonnegative(gtr in arb_gtr_nuc(), t in arb_branch_len()) {
      let p = gtr.expQt(t);
      prop_assert_array_nonneg!(p, epsilon = 1e-14);
    }

    #[test]
    fn test_prop_gtr_expqt_bounded(gtr in arb_gtr_nuc(), t in arb_branch_len()) {
      let p = gtr.expQt(t);
      prop_assert_array_upper_bounded!(p, bound = 1.0, epsilon = 1e-14);
    }

    #[test]
    fn test_prop_gtr_expqt_zero_is_identity(gtr in arb_gtr_nuc()) {
      let p = gtr.expQt(0.0);
      prop_assert_array_abs_diff_eq!(p, Array2::eye(4), epsilon = 1e-10);
    }

    #[test]
    fn test_prop_gtr_expqt_equilibrium_limit(gtr in arb_gtr_nuc()) {
      let t = 1000.0 / gtr.mu.max(0.001);
      let p = gtr.expQt(t);
      let expected = Array2::from_shape_fn((4, 4), |(i, _j)| gtr.pi[i]);
      prop_assert_array_abs_diff_eq!(p, expected, epsilon = 1e-6);
    }

    #[test]
    fn test_prop_gtr_expqt_semigroup(gtr in arb_gtr_nuc(), s in 0.001_f64..1.0, t in 0.001_f64..1.0) {
      let p_s = gtr.expQt(s);
      let p_t = gtr.expQt(t);
      let p_st = gtr.expQt(s + t);
      let p_s_times_p_t = p_s.dot(&p_t);
      prop_assert_array_abs_diff_eq!(p_st, p_s_times_p_t, epsilon = 1e-10);
    }

    #[test]
    fn test_prop_gtr_expqt_stationary_preserved(gtr in arb_gtr_nuc(), t in arb_branch_len()) {
      let p = gtr.expQt(t);
      let pi_evolved = p.dot(&gtr.pi);
      prop_assert_array_abs_diff_eq!(pi_evolved, gtr.pi, epsilon = 1e-10);
    }

    #[test]
    fn test_prop_gtr_expqt_evolve_transpose_of_propagate(
      gtr in arb_gtr_nuc(),
      profile in arb_profile_nuc(5),
      t in arb_branch_len()
    ) {
      let p = gtr.expQt(t);
      let p_t = p.t();

      let propagated = profile.dot(&p);
      let evolved = profile.dot(&p_t);

      let propagated_gtr = gtr.propagate_profile(&profile, t, false);
      let evolved_gtr = gtr.evolve(&profile, t, false);

      prop_assert_array_abs_diff_eq!(propagated_gtr, propagated, epsilon = 1e-10);

      prop_assert_array_abs_diff_eq!(evolved_gtr, evolved, epsilon = 1e-10);

      prop_assert_rows_sum_to(&evolved_gtr, 1.0, 1e-10)?;
    }
  }
}
