#[cfg(test)]
mod tests {
  use crate::gtr::__tests__::generators::tests::generators::{arb_branch_len, arb_gtr_nuc, arb_profile_nuc};
  use crate::gtr::__tests__::prop_support::prop_assert_rows_sum_to;
  use proptest::prelude::*;
  use treetime_utils::{prop_assert_array_finite, prop_assert_array_nonneg};

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    #[test]
    fn test_prop_gtr_numerical_expqt_finite(gtr in arb_gtr_nuc(), t in arb_branch_len()) {
      let p = gtr.expQt(t);
      prop_assert_array_finite!(p);
    }

    #[test]
    fn test_prop_gtr_numerical_propagate_profile_valid_output(
      gtr in arb_gtr_nuc(),
      profile in arb_profile_nuc(5),
      t in arb_branch_len()
    ) {
      let propagated = gtr.propagate_profile(&profile, t, false);

      prop_assert_array_finite!(propagated);
      prop_assert_array_nonneg!(propagated, epsilon = 1e-14);
    }

    #[test]
    fn test_prop_gtr_numerical_evolve_preserves_probability(
      gtr in arb_gtr_nuc(),
      profile in arb_profile_nuc(5),
      t in arb_branch_len()
    ) {
      let evolved = gtr.evolve(&profile, t, false);

      prop_assert_rows_sum_to(&evolved, 1.0, 1e-10)?;
      prop_assert_array_nonneg!(evolved, epsilon = 1e-14);
    }
  }
}
