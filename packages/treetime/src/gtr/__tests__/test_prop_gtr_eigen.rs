#[cfg(test)]
mod tests {
  use crate::gtr::__tests__::generators::tests::generators::arb_gtr_nuc;
  use ndarray::Array2;
  use proptest::prelude::*;
  use treetime_utils::{prop_assert_array_abs_diff_eq, prop_assert_array_upper_bounded};

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    #[test]
    fn test_prop_gtr_eigen_eigvals_nonpositive(gtr in arb_gtr_nuc()) {
      prop_assert_array_upper_bounded!(gtr.eigvals, bound = 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_prop_gtr_eigen_eigvals_one_zero_eigenvalue(gtr in arb_gtr_nuc()) {
      let zero_count = gtr.eigvals.iter().filter(|&&e| e.abs() < 1e-10).count();
      prop_assert!(
        zero_count == 1,
        "Expected exactly 1 zero eigenvalue, found {zero_count}: {:?}",
        gtr.eigvals
      );
    }

    #[test]
    fn test_prop_gtr_eigen_eigendecomposition_v_times_v_inv_is_identity(gtr in arb_gtr_nuc()) {
      let product = gtr.v.dot(&gtr.v_inv);
      prop_assert_array_abs_diff_eq!(product, Array2::eye(4), epsilon = 1e-10);
    }
  }
}
