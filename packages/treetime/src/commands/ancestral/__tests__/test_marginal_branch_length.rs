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
