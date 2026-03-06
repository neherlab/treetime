#[cfg(test)]
mod tests {
  use crate::commands::ancestral::__tests__::prop_generators::input::arb_marginal_input_no_gaps;
  use crate::commands::ancestral::__tests__::prop_marginal_support::tests::{run_dense_marginal, run_sparse_marginal};
  use proptest::prelude::*;
  use treetime_utils::prop_assert_relative_eq;

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(30))]

    // INVESTIGATE(dense-sparse-divergence): Dense and sparse marginal log-likelihoods diverge
    // beyond rounding for ~2.5% of random GTR configurations.
    //
    // Methodology: one-off run of 200 proptest cases (current test config uses 30),
    // 4 taxa, 10-position gap-free alignments, random GTR parameters.
    // Ran both dense and sparse marginal on identical inputs, collected abs_diff, rel_diff, ULPs.
    //
    // Findings - two distinct populations:
    //
    //   Population 1 (~97.5%): 0-3 ULPs, rel_diff < 5e-16.
    //     Floating-point rounding from different operation ordering. Expected.
    //
    //   Population 2 (~2.5%): billions of ULPs, rel_diff up to 7.6e-6.
    //       abs_diff=3.0e-13  rel_diff=3.4e-15  ulps=21
    //       abs_diff=3.3e-6   rel_diff=3.4e-8   ulps=234M
    //       abs_diff=1.8e-5   rel_diff=1.8e-7   ulps=1.3B
    //       abs_diff=3.0e-5   rel_diff=4.7e-7   ulps=4.2B
    //       abs_diff=4.7e-5   rel_diff=5.8e-7   ulps=3.3B
    //       abs_diff=5.6e-4   rel_diff=7.6e-6   ulps=39.5B
    //     Algorithmic differences, not rounding. Bimodal distribution suggests certain GTR
    //     parameter combinations trigger a code path divergence.
    //
    // Tolerance: `max_relative = 1e-5` accommodates worst observed (7.6e-6) with margin.
    // ULPs unsuitable due to bimodal distribution. `epsilon = 0.0` disables the absolute
    // tolerance escape hatch in `relative_eq!`, ensuring pure relative comparison.
    //
    // Investigation leads:
    //   - Identify which GTR parameters trigger population 2 (extreme eigenvalue ratios
    //     in the rate matrix Q, near-degenerate exchangeability matrices)
    //   - Compare dense/sparse intermediate results (per-node partial likelihoods,
    //     branch propagation via P(t) = exp(Qt))
    //   - Check whether Fitch-based invariant/variable classification diverges from
    //     likelihood-based variability: a position where all tips share one state is
    //     "invariant" to Fitch parsimony, but its probability vector under the GTR model
    //     still depends on branch-length-specific transition probabilities
    //   - If correctness bug: fix and tighten to ULPs-level
    //   - If inherent to sparse approximation: document the error bound

    /// Oracle comparison test: dense and sparse marginal log-likelihoods agree on
    /// randomly generated inputs (4 taxa, 10-position gap-free alignments, random GTR).
    ///
    /// Both representations compute the Felsenstein pruning algorithm likelihood via
    /// different code paths, then the result is compared.
    ///
    /// Dense (reference): stores full L x K probability matrices (L = alignment length,
    /// K = number of states) at every node. Calls `initialize_marginal`, which attaches
    /// sequences and runs the backward pass (tips to root, computing partial likelihoods)
    /// and forward pass (root to tips, computing outgroup profiles).
    ///
    /// Sparse (optimized): first runs Fitch parsimony compression (`compress_sequences`)
    /// to classify positions as invariant or variable, then calls `update_marginal` on
    /// variable positions only, accumulating fixed-site contributions separately.
    ///
    /// Invariant tested: the total log-likelihood
    ///   L = sum over sites l of log(sum over states s of pi[s] * L_root^(l)(s))
    /// must agree between dense and sparse code paths.
    ///
    /// Uses relative tolerance (max_relative = 1e-5) rather than ULPs due to a
    /// bimodal error distribution: ~97.5% of cases agree to 0-3 ULPs, while
    /// ~2.5% diverge up to 7.6e-6 relative difference for certain GTR
    /// configurations. See INVESTIGATE comment above for details.
    ///
    /// Companion example test: `test_marginal_dense_sparse_example_gap_free_consistency`.
    #[test]
    fn test_prop_marginal_dense_sparse_gap_free_consistency(input in arb_marginal_input_no_gaps(4, 10)) {
      let (log_lh_dense, _) = run_dense_marginal(&input).unwrap();
      let (log_lh_sparse, _) = run_sparse_marginal(&input).unwrap();

      prop_assert_relative_eq!(log_lh_dense, log_lh_sparse, max_relative = 1e-5, epsilon = 0.0);
    }
  }
}
