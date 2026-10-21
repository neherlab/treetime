#[cfg(test)]
mod tests {
  use crate::ancestral::__tests__::prop_generators::input::arb_marginal_input_no_gaps;
  use crate::ancestral::__tests__::prop_marginal_support::tests::{
    run_dense_marginal_for_input, run_sparse_marginal_for_input,
  };
  use proptest::prelude::*;
  use treetime_utils::prop_assert_relative_eq;

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(30))]

    // TODO: kb/issues/M-ancestral-dense-sparse-divergence.md
    #[test]
    #[ignore = "dense and sparse marginal log-likelihoods diverge for some GTR models: kb/issues/M-ancestral-dense-sparse-divergence.md"]
    fn test_prop_marginal_dense_sparse_gap_free_consistency(input in arb_marginal_input_no_gaps(4, 10)) {
      let (log_lh_dense, _) = run_dense_marginal_for_input(&input).unwrap();
      let (log_lh_sparse, _) = run_sparse_marginal_for_input(&input).unwrap();

      prop_assert_relative_eq!(log_lh_dense, log_lh_sparse, max_relative = 1e-5, epsilon = 0.0);
    }
  }
}
