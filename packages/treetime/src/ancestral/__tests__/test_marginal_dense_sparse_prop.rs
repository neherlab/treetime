#[cfg(test)]
mod tests {
  use crate::ancestral::__tests__::prop_generators::input::arb_marginal_input_no_gaps;
  use crate::ancestral::__tests__::prop_marginal_support::tests::{run_dense_marginal, run_sparse_marginal};
  use proptest::prelude::*;
  use treetime_utils::prop_assert_relative_eq;

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(30))]

    #[test]
    #[ignore = "Investigate dense-sparse marginal log-likelihood divergence on certain GTR configs (see INVESTIGATE comment)"]
    fn test_prop_marginal_dense_sparse_gap_free_consistency(input in arb_marginal_input_no_gaps(4, 10)) {
      let (log_lh_dense, _) = run_dense_marginal(&input).unwrap();
      let (log_lh_sparse, _) = run_sparse_marginal(&input).unwrap();

      prop_assert_relative_eq!(log_lh_dense, log_lh_sparse, max_relative = 1e-5, epsilon = 0.0);
    }
  }
}
