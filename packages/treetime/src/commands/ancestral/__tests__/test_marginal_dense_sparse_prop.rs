#[cfg(test)]
mod tests {
  use crate::commands::ancestral::__tests__::prop_generators::input::arb_marginal_input_no_gaps;
  use crate::commands::ancestral::__tests__::prop_marginal_support::{run_dense_marginal, run_sparse_marginal};
  use proptest::prelude::*;

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(30))]

    /// Companion example test: `test_marginal_dense_sparse_example_gap_free_consistency`.
    #[test]
    fn test_prop_marginal_dense_sparse_gap_free_consistency(input in arb_marginal_input_no_gaps(4, 10)) {
      let (log_lh_dense, _) = run_dense_marginal(&input).unwrap();
      let (log_lh_sparse, _) = run_sparse_marginal(&input).unwrap();

      let abs_diff = (log_lh_dense - log_lh_sparse).abs();
      let scale = log_lh_dense.abs().max(log_lh_sparse.abs()).max(1.0);
      let rel_diff = abs_diff / scale;
      prop_assert!(
        rel_diff < 1e-6,
        "Dense/sparse mismatch: dense={log_lh_dense}, sparse={log_lh_sparse}, rel_diff={rel_diff}"
      );
    }
  }
}
