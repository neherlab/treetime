#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::gtr::gtr::{GTR, GTRParams};
  use crate::pretty_assert_ulps_eq;
  use crate::test_utils::{run_dense_marginal_with_newick, run_sparse_marginal_with_newick};
  use eyre::Report;
  use itertools::Itertools;
  use ndarray::prelude::*;

  /// Small epsilon for monotonicity comparisons.
  /// Allows for floating-point noise while catching algorithmic errors.
  const MONOTONICITY_EPSILON: f64 = 1e-10;

  /// Epsilon for equilibrium convergence tests.
  const EQUILIBRIUM_EPSILON: f64 = 1e-6;

  // ============================================================================
  // T7: Likelihood monotonicity tests
  // ============================================================================

  /// Likelihood increases monotonically with branch length for mismatched sequences.
  ///
  /// Tree: (A:t, B:t)root; leaf A='A', leaf B='T'. Branch length t varies from 0.1 to 2.0.
  ///
  /// At t=0, the transition matrix is the identity, so P(T|A,0)=0 and the mismatch is
  /// impossible (likelihood near zero). As t increases, off-diagonal elements of expQt
  /// grow, making the required A->T or T->A transition more probable.
  ///
  /// At t->infinity, expQt rows converge to pi, so L -> pi[A] * pi[T] = 0.0625 under JC69.
  /// The test verifies monotonic increase and convergence toward equilibrium.
  #[test]
  fn test_likelihood_monotonic_increase_mismatched_sequences() -> Result<(), Report> {
    // For mismatched sequences (A vs T), likelihood increases monotonically
    // from t=0 to t=infinity:
    // - At t=0: very low (no time for transition, mismatch impossible)
    // - As t increases: more likely that transitions occurred
    // - At t->infinity: converges to equilibrium pi[A] * pi[T] = 0.0625
    //
    // Unlike matched sequences, mismatched sequences BENEFIT from longer branches
    // because the probability of the required transition increases.
    let gtr = jc69(JC69Params::default())?;

    let branch_lengths: Vec<f64> = (1..=20).map(|i| i as f64 * 0.1).collect_vec();
    let aln = ">A\nA\n>B\nT\n";

    let mut prev_log_lh = f64::NEG_INFINITY;
    for t in &branch_lengths {
      let newick = format!("(A:{t},B:{t})root;");
      let log_lh = run_dense_marginal_with_newick(&newick, aln, gtr.clone())?;

      // Verify monotonic increase (likelihood increases as branch length increases)
      assert!(
        log_lh >= prev_log_lh - MONOTONICITY_EPSILON,
        "Likelihood should increase monotonically for mismatched sequences. \
         At t={t}: log_lh={log_lh}, prev={prev_log_lh}"
      );
      prev_log_lh = log_lh;
    }

    // Verify it converges towards equilibrium
    let equilibrium_log_lh = (0.25 * 0.25_f64).ln();
    let final_log_lh = prev_log_lh;
    assert!(
      (final_log_lh - equilibrium_log_lh).abs() < 0.1,
      "At t=2.0, should be close to equilibrium. \
       actual={final_log_lh}, equilibrium={equilibrium_log_lh}"
    );

    Ok(())
  }

  /// Likelihood decreases monotonically with branch length for matched sequences.
  ///
  /// Tree: (A:t, B:t)root; both leaves observe 'A'. Branch length t varies from 0.1 to 2.0.
  ///
  /// At t=0, expQt is the identity matrix, so P(A|A,0)=1 and the matched observation
  /// has maximum probability. As t increases, off-diagonal elements grow, reducing the
  /// probability that both leaves independently retain the ancestral state.
  ///
  /// This is the dual of the mismatched test: matched sequences are best explained
  /// by short branches (few substitutions), while longer branches make the match
  /// increasingly unlikely.
  #[test]
  fn test_likelihood_maximized_near_zero_for_matched_sequences() -> Result<(), Report> {
    // For matched sequences (A and A), likelihood should be maximized near t=0
    // and decrease monotonically as branch length increases.
    let gtr = jc69(JC69Params::default())?;

    let branch_lengths: Vec<f64> = (1..=20).map(|i| i as f64 * 0.1).collect_vec();

    let mut prev_log_lh = f64::INFINITY;
    for t in &branch_lengths {
      let newick = format!("(A:{t},B:{t})root;");
      let aln = ">A\nA\n>B\nA\n";
      let log_lh = run_dense_marginal_with_newick(&newick, aln, gtr.clone())?;

      // For matched sequences, likelihood monotonically decreases as t increases
      assert!(
        log_lh <= prev_log_lh + MONOTONICITY_EPSILON,
        "Likelihood should decrease as branch length increases for matched sequences. \
         At t={t}: log_lh={log_lh}, prev={prev_log_lh}"
      );
      prev_log_lh = log_lh;
    }
    Ok(())
  }

  /// Numerical stability: log-likelihood remains finite and non-positive across 4 orders
  /// of magnitude in branch length (0.001 to 10.0).
  ///
  /// Tree: (A:t, B:t)root; sequences "ACGT" vs "TGCA" (all mismatched).
  ///
  /// At extreme branch lengths, matrix exponentiation (expQt) can produce values near
  /// machine epsilon (short t, off-diagonal) or near equilibrium (long t). This test
  /// verifies that no intermediate computation produces NaN or infinity, and that
  /// log-likelihood remains in the valid range (-inf, 0].
  #[test]
  fn test_likelihood_finite_across_branch_length_range() -> Result<(), Report> {
    // Verify no NaN/Inf across a wide range of branch lengths
    let gtr = jc69(JC69Params::default())?;

    let branch_lengths = [0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0];

    for t in branch_lengths {
      let newick = format!("(A:{t},B:{t})root;");
      let aln = ">A\nACGT\n>B\nTGCA\n";
      let log_lh = run_dense_marginal_with_newick(&newick, aln, gtr.clone())?;

      assert!(
        log_lh.is_finite(),
        "Log-likelihood should be finite at t={t}, got {log_lh}"
      );
      assert!(
        log_lh <= 0.0,
        "Log-likelihood should be non-positive at t={t}, got {log_lh}"
      );
    }
    Ok(())
  }

  /// Monotonic decrease of likelihood for identical sequences on a three-taxon tree.
  ///
  /// Tree: ((A:t, B:t)AB:t, C:t)root; all branches have the same length t.
  /// All three leaves observe "AAAA".
  ///
  /// With identical sequences, the maximum-likelihood branch length is zero. As all
  /// branches grow uniformly, the probability of preserving state 'A' along every path
  /// decreases monotonically. This tests monotonicity on a non-trivial tree topology
  /// where messages must propagate through the internal node AB before reaching the root.
  #[test]
  fn test_likelihood_monotonicity_three_taxon_tree() -> Result<(), Report> {
    // Test monotonicity on a more complex tree structure
    // Tree: ((A:t,B:t)AB:t,C:t)root; with all same branch length
    let gtr = jc69(JC69Params::default())?;

    // All sequences same - likelihood should decrease with branch length
    let branch_lengths: Vec<f64> = (1..=15).map(|i| i as f64 * 0.1).collect_vec();

    let mut prev_log_lh = f64::INFINITY;
    for t in &branch_lengths {
      let newick = format!("((A:{t},B:{t})AB:{t},C:{t})root;");
      let aln = ">A\nAAAA\n>B\nAAAA\n>C\nAAAA\n";
      let log_lh = run_dense_marginal_with_newick(&newick, aln, gtr.clone())?;

      assert!(
        log_lh <= prev_log_lh + MONOTONICITY_EPSILON,
        "Likelihood should decrease for identical sequences as t increases. \
         At t={t}: log_lh={log_lh}, prev={prev_log_lh}"
      );
      prev_log_lh = log_lh;
    }
    Ok(())
  }

  /// Monotonic decrease of likelihood for identical sequences using sparse representation.
  ///
  /// Same invariant as `test_likelihood_maximized_near_zero_for_matched_sequences`, but
  /// using the sparse partition path (`run_sparse_marginal_with_newick`) instead of dense.
  /// Since all positions are identical across both leaves, the sparse representation compresses
  /// the alignment to zero variable positions, exercising the invariant-position likelihood
  /// accumulation code path.
  #[test]
  fn test_likelihood_monotonicity_sparse_partition() -> Result<(), Report> {
    // Verify same monotonicity behavior with sparse partition
    let gtr = jc69(JC69Params::default())?;

    let branch_lengths: Vec<f64> = (1..=15).map(|i| i as f64 * 0.1).collect_vec();

    let mut prev_log_lh = f64::INFINITY;
    for t in &branch_lengths {
      let newick = format!("(A:{t},B:{t})root;");
      let aln = ">A\nGGGG\n>B\nGGGG\n";
      let log_lh = run_sparse_marginal_with_newick(&newick, aln, gtr.clone())?;

      assert!(
        log_lh <= prev_log_lh + MONOTONICITY_EPSILON,
        "Sparse: Likelihood should decrease for identical sequences. \
         At t={t}: log_lh={log_lh}, prev={prev_log_lh}"
      );
      prev_log_lh = log_lh;
    }
    Ok(())
  }

  // ============================================================================
  // T8: Equilibrium convergence tests
  // ============================================================================

  /// At very long branch lengths, likelihood converges to the product of equilibrium
  /// frequencies (dense partition).
  ///
  /// Tree: (A:100, B:100)root; leaf A='A', leaf B='T'.
  ///
  /// As t -> infinity, each row of expQt converges to pi (the equilibrium distribution).
  /// The transition probability becomes independent of the ancestral state:
  ///   P(i|s, t->inf) -> pi[i]
  ///
  /// Therefore:
  ///   L -> sum_s pi[s] * pi[obs_A] * pi[obs_B] = pi[obs_A] * pi[obs_B]
  ///
  /// Under JC69: L -> 0.25 * 0.25 = 0.0625, so ln(L) = ln(0.0625) ~ -2.77.
  #[test]
  fn test_equilibrium_convergence_dense() -> Result<(), Report> {
    // As branch length -> infinity, the transition matrix converges to
    // a matrix where each row equals pi (equilibrium frequencies).
    // For a two-taxon tree with very long branches, the likelihood should
    // approach: sum_s pi[s] * pi[obs_A] * pi[obs_B] = pi[obs_A] * pi[obs_B]
    let gtr = jc69(JC69Params::default())?;
    let t = 100.0; // Very long branch

    // For JC69, pi = [0.25, 0.25, 0.25, 0.25]
    // At equilibrium: L ≈ pi[A] * pi[T] = 0.25 * 0.25 = 0.0625
    // log(0.0625) ≈ -2.77
    let expected_equilibrium_log_lh = (0.25 * 0.25_f64).ln();

    let newick = format!("(A:{t},B:{t})root;");
    let aln = ">A\nA\n>B\nT\n";
    let actual_log_lh = run_dense_marginal_with_newick(&newick, aln, gtr)?;

    pretty_assert_ulps_eq!(
      expected_equilibrium_log_lh,
      actual_log_lh,
      epsilon = EQUILIBRIUM_EPSILON
    );
    Ok(())
  }

  /// Equilibrium convergence test using sparse partition. Same invariant as
  /// `test_equilibrium_convergence_dense`: at t=100, L -> pi[A] * pi[T] = 0.0625.
  ///
  /// Verifies that the sparse representation's handling of variable positions produces
  /// the same equilibrium limit as the dense path.
  #[test]
  fn test_equilibrium_convergence_sparse() -> Result<(), Report> {
    // Same equilibrium test with sparse partition
    let gtr = jc69(JC69Params::default())?;
    let t = 100.0;

    let expected_equilibrium_log_lh = (0.25 * 0.25_f64).ln();

    let newick = format!("(A:{t},B:{t})root;");
    let aln = ">A\nA\n>B\nT\n";
    let actual_log_lh = run_sparse_marginal_with_newick(&newick, aln, gtr)?;

    pretty_assert_ulps_eq!(
      expected_equilibrium_log_lh,
      actual_log_lh,
      epsilon = EQUILIBRIUM_EPSILON
    );
    Ok(())
  }

  /// Equilibrium convergence with non-uniform equilibrium frequencies (GTR model).
  ///
  /// Uses pi = [0.4, 0.1, 0.2, 0.3]; leaf A='A', leaf B='G'.
  /// At t=100: L -> pi[A] * pi[G] = 0.4 * 0.2 = 0.08.
  ///
  /// Non-uniform pi changes the equilibrium limit. This verifies that the implementation
  /// correctly propagates the model's equilibrium frequencies into the long-branch limit,
  /// not just the JC69 uniform case.
  #[test]
  fn test_equilibrium_convergence_nonuniform_pi() -> Result<(), Report> {
    // Test equilibrium convergence with non-uniform equilibrium frequencies
    let alphabet = Alphabet::new(AlphabetName::Nuc, true)?;
    let gtr = GTR::new(GTRParams {
      alphabet,
      mu: 1.0,
      W: None,
      pi: array![0.4, 0.1, 0.2, 0.3], // A=0.4, C=0.1, G=0.2, T=0.3
    })?;

    let t = 100.0;

    // At equilibrium with A at one leaf, G at other:
    // L ≈ pi[A] * pi[G] = 0.4 * 0.2 = 0.08
    let expected_equilibrium_log_lh = (0.4 * 0.2_f64).ln();

    let newick = format!("(A:{t},B:{t})root;");
    let aln = ">A\nA\n>B\nG\n";
    let actual_log_lh = run_dense_marginal_with_newick(&newick, aln, gtr)?;

    pretty_assert_ulps_eq!(
      expected_equilibrium_log_lh,
      actual_log_lh,
      epsilon = EQUILIBRIUM_EPSILON
    );
    Ok(())
  }

  /// Multi-position equilibrium convergence under JC69.
  ///
  /// Sequences: A="ACG", B="TCA" (3 positions, all mismatched). At t=100, each position
  /// independently converges to pi[obs_A_i] * pi[obs_B_i] = 0.25 * 0.25 = 0.0625.
  ///
  /// Total: ln(L) = 3 * ln(0.0625) ~ -8.32.
  ///
  /// Verifies that per-position equilibrium limits compose correctly via log-sum when
  /// the alignment has multiple independent columns.
  #[test]
  fn test_equilibrium_convergence_multiple_positions() -> Result<(), Report> {
    // Multi-position equilibrium test
    let gtr = jc69(JC69Params::default())?;
    let t = 100.0;

    // Sequence: ACG vs TCA
    // At equilibrium: product of pi[obs_A_i] * pi[obs_B_i] for each position
    // = (0.25 * 0.25) * (0.25 * 0.25) * (0.25 * 0.25) = 0.0625^3
    let expected_equilibrium_log_lh = 3.0 * (0.25 * 0.25_f64).ln();

    let newick = format!("(A:{t},B:{t})root;");
    let aln = ">A\nACG\n>B\nTCA\n";
    let actual_log_lh = run_dense_marginal_with_newick(&newick, aln, gtr)?;

    pretty_assert_ulps_eq!(
      expected_equilibrium_log_lh,
      actual_log_lh,
      epsilon = EQUILIBRIUM_EPSILON
    );
    Ok(())
  }

  /// Equilibrium convergence on a 4-leaf star tree under JC69.
  ///
  /// Tree: (A:100, B:100, C:100, D:100)root; leaves observe A, C, G, T.
  /// At equilibrium: L -> pi[A] * pi[C] * pi[G] * pi[T] = 0.25^4 = 1/256.
  ///
  /// With 4 leaves, the product of 4 equilibrium terms yields a very small likelihood.
  /// Verifies that the n-ary root message-passing correctly accumulates the equilibrium
  /// limit across all children.
  #[test]
  fn test_equilibrium_star_tree() -> Result<(), Report> {
    // Star tree with long branches should also converge to equilibrium
    let gtr = jc69(JC69Params::default())?;
    let t = 100.0;

    // 4 leaves with states A, C, G, T
    // At equilibrium: pi[A] * pi[C] * pi[G] * pi[T] = 0.25^4
    let expected_equilibrium_log_lh = (0.25_f64.powi(4)).ln();

    let newick = format!("(A:{t},B:{t},C:{t},D:{t})root;");
    let aln = ">A\nA\n>B\nC\n>C\nG\n>D\nT\n";
    let actual_log_lh = run_dense_marginal_with_newick(&newick, aln, gtr)?;

    pretty_assert_ulps_eq!(
      expected_equilibrium_log_lh,
      actual_log_lh,
      epsilon = EQUILIBRIUM_EPSILON
    );
    Ok(())
  }

  /// Numerical stability at extremely short branch lengths (t=1e-8).
  ///
  /// Tree: (A:1e-8, B:1e-8)root; both leaves observe 'A'.
  ///
  /// At t -> 0, expQt -> I (identity matrix), so P(A|s,t) -> delta(A,s). The likelihood
  /// simplifies to L -> pi[A] = 0.25, and ln(L) -> ln(0.25) ~ -1.386.
  ///
  /// Tests that the matrix exponentiation and likelihood computation remain finite and
  /// accurate when branch lengths approach machine epsilon. Underflow in expQt off-diagonal
  /// elements or loss of precision in the diagonal could produce NaN or incorrect values.
  #[test]
  fn test_branch_length_sensitivity_near_zero() -> Result<(), Report> {
    // Near t=0 with matched sequences, the likelihood should be close to
    // the analytical formula from test_marginal_analytical.rs.
    // This test verifies behavior at extreme short branch lengths.
    let gtr = jc69(JC69Params::default())?;
    let t = 1e-8;

    let newick = format!("(A:{t},B:{t})root;");
    let aln = ">A\nA\n>B\nA\n";
    let log_lh = run_dense_marginal_with_newick(&newick, aln, gtr)?;

    // At very short branches with matched sequences:
    // - Log-likelihood should be finite and negative
    // - Should be close to ln(pi[A]) = ln(0.25) ≈ -1.386 since expQt ≈ I
    assert!(log_lh.is_finite(), "Log-likelihood should be finite at t={t}");
    assert!(log_lh <= 0.0, "Log-likelihood should be non-positive at t={t}");

    // Should be close to ln(0.25) with some tolerance for implementation details
    let approx_expected = 0.25_f64.ln(); // -1.386...
    assert!(
      (log_lh - approx_expected).abs() < 0.01,
      "At very short branches, log_lh should be close to ln(pi). \
       Expected ≈{approx_expected}, got {log_lh}"
    );

    Ok(())
  }

  /// Dense and sparse partitions produce identical log-likelihoods across branch lengths.
  ///
  /// Tree: (A:t, B:t)root; sequences "ACGTACGT" vs "TGCATGCA" (all mismatched).
  /// Branch lengths: 0.01, 0.1, 0.5, 1.0, 5.0, 20.0.
  ///
  /// The dense path stores full probability vectors at every position. The sparse path
  /// compresses invariant positions and only computes marginals at variable sites. Both
  /// paths implement the same Felsenstein pruning algorithm and must produce identical
  /// results (within floating-point tolerance of 1e-10).
  ///
  /// Testing across a wide range of branch lengths exercises both paths under different
  /// numerical regimes: near-identity matrices (small t), moderate substitution rates
  /// (medium t), and near-equilibrium (large t).
  #[test]
  fn test_dense_sparse_consistency_across_branch_lengths() -> Result<(), Report> {
    // Dense and sparse partitions should give same likelihoods across branch lengths
    let gtr = jc69(JC69Params::default())?;

    let branch_lengths = [0.01, 0.1, 0.5, 1.0, 5.0, 20.0];

    for t in branch_lengths {
      let newick = format!("(A:{t},B:{t})root;");
      let aln = ">A\nACGTACGT\n>B\nTGCATGCA\n";

      let dense_log_lh = run_dense_marginal_with_newick(&newick, aln, gtr.clone())?;
      let sparse_log_lh = run_sparse_marginal_with_newick(&newick, aln, gtr.clone())?;

      pretty_assert_ulps_eq!(
        dense_log_lh,
        sparse_log_lh,
        epsilon = 1e-10,
        "Dense and sparse should match at t={t}"
      );
    }
    Ok(())
  }
}
