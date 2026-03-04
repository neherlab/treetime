#[cfg(test)]
mod tests {
  //! Property-based numerical stability tests for GTR models.
  //!
  //! Verifies no NaN/Inf values and valid probability outputs across
  //! random GTR models and branch lengths.

  use crate::gtr::__tests__::generators::tests::generators::{arb_branch_len, arb_gtr_nuc, arb_profile_nuc};
  use crate::gtr::__tests__::prop_support::prop_assert_rows_sum_to;
  use proptest::prelude::*;

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    /// NaN typically arises from:
    /// - 0/0 in normalization
    /// - sqrt of negative number
    /// - inf - inf in eigendecomposition
    /// - log(0) or log(negative) if using log-space
    ///
    /// Test across random valid GTR models and branch lengths to ensure
    /// no input combination produces NaN.
    #[test]
    fn test_prop_gtr_numerical_no_nan_in_expqt(gtr in arb_gtr_nuc(), t in arb_branch_len()) {
      let p = gtr.expQt(t);
      for i in 0..4 {
        for j in 0..4 {
          prop_assert!(
            !p[[i, j]].is_nan(),
            "P(t={t})[{i},{j}] is NaN"
          );
        }
      }
    }

    /// Inf typically arises from:
    /// - Overflow in exp(large positive number) - shouldn't happen with valid Q
    /// - Division by very small number in V^{-1}
    /// - Accumulation of large values in matrix multiplication
    ///
    /// Test across random valid GTR models to ensure no Inf values.
    #[test]
    fn test_prop_gtr_numerical_no_inf_in_expqt(gtr in arb_gtr_nuc(), t in arb_branch_len()) {
      let p = gtr.expQt(t);
      for i in 0..4 {
        for j in 0..4 {
          prop_assert!(
            !p[[i, j]].is_infinite(),
            "P(t={t})[{i},{j}] is Inf"
          );
        }
      }
    }

    /// propagate_profile computes p @ P(t).
    /// In transposed convention, P is column-stochastic, so p @ P doesn't
    /// preserve row sums (unlike p @ P.T which does for evolve).
    ///
    /// However, the output should still be valid (non-negative, no NaN/Inf).
    #[test]
    fn test_prop_gtr_numerical_propagate_profile_valid_output(
      gtr in arb_gtr_nuc(),
      profile in arb_profile_nuc(5),
      t in arb_branch_len()
    ) {
      let propagated = gtr.propagate_profile(&profile, t, false);

      // Check no NaN or Inf
      for i in 0..5 {
        for j in 0..4 {
          prop_assert!(
            !propagated[[i, j]].is_nan(),
            "propagate_profile[{i},{j}] is NaN at t={t}"
          );
          prop_assert!(
            !propagated[[i, j]].is_infinite(),
            "propagate_profile[{i},{j}] is Inf at t={t}"
          );
          prop_assert!(
            propagated[[i, j]] >= -1e-14,
            "propagate_profile[{i},{j}] = {} is negative at t={t}",
            propagated[[i, j]]
          );
        }
      }
    }

    /// evolve computes p @ P(t).T.
    /// In transposed convention, P.T is row-stochastic, so evolve DOES
    /// preserve row sums: output rows should sum to 1.
    #[test]
    fn test_prop_gtr_numerical_evolve_preserves_probability(
      gtr in arb_gtr_nuc(),
      profile in arb_profile_nuc(5),
      t in arb_branch_len()
    ) {
      let evolved = gtr.evolve(&profile, t, false);

      prop_assert_rows_sum_to(&evolved, 1.0, 1e-10)?;
      for ((i, j), &val) in evolved.indexed_iter() {
        prop_assert!(val >= -1e-14, "evolve[{i},{j}] = {val} is negative");
      }
    }
  }
}
