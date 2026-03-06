#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::commands::ancestral::fitch::{compress_sequences, get_common_length};
  use crate::commands::ancestral::marginal::{initialize_marginal, update_marginal};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::gtr::gtr::{GTR, GTRParams};
  use crate::representation::partition::marginal_dense::PartitionMarginalDense;
  use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
  use crate::representation::payload::ancestral::GraphAncestral;
  use crate::representation::payload::dense::DenseSeqDis;
  use crate::representation::payload::sparse::MarginalSparseSeqDistribution;
  use approx::assert_ulps_eq;
  use eyre::Report;
  use maplit::btreemap;
  use ndarray::prelude::*;
  use parking_lot::RwLock;
  use std::sync::{Arc, LazyLock};
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;

  static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

  /// Assert that a dense marginal profile is numerically stable.
  ///
  /// Checks four conditions:
  /// - `log_lh` is finite.
  /// - Row sum equals 1.0 within `max_ulps` units in last place (at every position).
  /// - All values are finite (no NaN or Inf from overflow/underflow).
  /// - All values are non-negative (probabilities in [0, 1]).
  ///
  /// The ULP-based tolerance scales with floating-point magnitude, making it
  /// appropriate for detecting numerical drift across diverse parameter regimes.
  fn assert_dense_profile_stable(profile: &DenseSeqDis, max_ulps: u32) {
    assert!(
      profile.log_lh.is_finite(),
      "Profile log_lh is not finite: {}",
      profile.log_lh
    );

    for (pos, row) in profile.dis.outer_iter().enumerate() {
      let sum: f64 = row.sum();
      assert_ulps_eq!(sum, 1.0, max_ulps = max_ulps);
      for (idx, &val) in row.iter().enumerate() {
        assert!(
          val.is_finite(),
          "Position {pos}, index {idx} has non-finite value: {val}"
        );
        assert!(val >= -1e-15, "Position {pos}, index {idx} has negative value: {val}");
      }
    }
  }

  /// Assert that a sparse marginal profile is numerically stable.
  ///
  /// Checks the same conditions as `assert_dense_profile_stable` but adapted for
  /// the sparse representation:
  /// - `log_lh` is finite.
  /// - Every variable-position distribution sums to 1.0 within `max_ulps`.
  /// - Every fixed-character distribution sums to 1.0 within `max_ulps`.
  /// - All values are finite and non-negative.
  fn assert_sparse_profile_stable(profile: &MarginalSparseSeqDistribution, max_ulps: u32) {
    assert!(
      profile.log_lh.is_finite(),
      "Profile log_lh is not finite: {}",
      profile.log_lh
    );

    for (pos, var_pos) in &profile.variable {
      let sum: f64 = var_pos.dis.sum();
      assert_ulps_eq!(sum, 1.0, max_ulps = max_ulps);
      for (idx, &val) in var_pos.dis.iter().enumerate() {
        assert!(
          val.is_finite(),
          "Variable position {pos}, index {idx} has non-finite value: {val}"
        );
        assert!(
          val >= -1e-15,
          "Variable position {pos}, index {idx} has negative value: {val}"
        );
      }
    }

    for (char_key, fixed_dis) in &profile.fixed {
      let sum: f64 = fixed_dis.sum();
      assert_ulps_eq!(sum, 1.0, max_ulps = max_ulps);
      for (idx, &val) in fixed_dis.iter().enumerate() {
        assert!(
          val.is_finite(),
          "Fixed distribution for char {char_key:?}, index {idx} has non-finite value: {val}"
        );
        assert!(
          val >= -1e-15,
          "Fixed distribution for char {char_key:?}, index {idx} has negative value: {val}"
        );
      }
    }
  }

  /// Run dense marginal reconstruction with a custom GTR model and return the
  /// log-likelihood and partition array.
  ///
  /// Parses the Newick tree and FASTA alignment, constructs a single dense partition
  /// with the given GTR model, and runs `initialize_marginal` (two-pass marginal
  /// ancestral reconstruction: backward pass computes partial likelihoods via
  /// Felsenstein's pruning algorithm, forward pass propagates root information
  /// back to compute marginal posteriors). Returns the computed log-likelihood
  /// and the partition for inspection.
  fn run_dense_marginal_with_partitions(
    newick: &str,
    aln_str: &str,
    gtr: GTR,
  ) -> Result<(f64, [Arc<RwLock<PartitionMarginalDense>>; 1]), Report> {
    let graph: GraphAncestral = nwk_read_str(newick)?;
    let aln = read_many_fasta_str(aln_str, &*NUC_ALPHABET)?;
    let alphabet = Alphabet::new(AlphabetName::Nuc, true)?;

    let partitions = [Arc::new(RwLock::new(PartitionMarginalDense {
      index: 0,
      gtr,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    let log_lh = initialize_marginal(&graph, &partitions, &aln)?;
    Ok((log_lh, partitions))
  }

  /// Run sparse marginal reconstruction with a custom GTR model and return the
  /// log-likelihood and partition array.
  ///
  /// Parses the Newick tree and FASTA alignment, constructs a single sparse partition,
  /// runs Fitch compression (`compress_sequences`) to identify variable vs fixed
  /// positions, then runs `update_marginal`. Returns the computed log-likelihood
  /// and the partition for inspection.
  fn run_sparse_marginal_with_partitions(
    newick: &str,
    aln_str: &str,
    gtr: GTR,
  ) -> Result<(f64, [Arc<RwLock<PartitionMarginalSparse>>; 1]), Report> {
    let graph: GraphAncestral = nwk_read_str(newick)?;
    let aln = read_many_fasta_str(aln_str, &*NUC_ALPHABET)?;
    let alphabet = Alphabet::new(AlphabetName::Nuc, true)?;

    let partitions = [Arc::new(RwLock::new(PartitionMarginalSparse {
      index: 0,
      gtr,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    compress_sequences(&graph, &partitions, &aln)?;
    let log_lh = update_marginal(&graph, &partitions)?;
    Ok((log_lh, partitions))
  }

  // ============================================================================
  // T4: Extreme branch length tests
  // ============================================================================

  /// Numerical stability with extremely short branches (t = 1e-10), dense representation.
  ///
  /// When t is near zero, the transition probability matrix approaches the identity:
  /// P(t) ~ I + Qt. The implementation computes P(t) = v * diag(exp(lambda_i * mu * t)) * v_inv,
  /// where v and v_inv incorporate D^{+-1/2} factors from the symmetrization of the rate
  /// matrix (see `eig_single_site()`). When exp(lambda_i * mu * t) ~ 1 for all eigenvalues,
  /// the off-diagonal elements of P(t) are O(mu*t) ~ 1e-10, testing whether the
  /// back-transform preserves normalization at near-machine-epsilon scale.
  ///
  /// Uses identical sequences (ACGT) on both leaves so the likelihood is dominated
  /// by the diagonal of P(t).
  #[test]
  fn test_extreme_short_branch_dense() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let newick = "(A:1e-10,B:1e-10)root;";
    let aln = ">A\nACGT\n>B\nACGT\n";

    let (log_lh, partitions) = run_dense_marginal_with_partitions(newick, aln, gtr)?;

    assert!(log_lh.is_finite(), "Log-likelihood is not finite: {log_lh}");
    assert!(log_lh <= 0.0, "Log-likelihood should be non-positive: {log_lh}");

    let partition = partitions[0].read_arc();
    for node_data in partition.nodes.values() {
      assert_dense_profile_stable(&node_data.profile, 8);
    }

    Ok(())
  }

  /// Numerical stability with extremely short branches (t = 1e-10), sparse representation.
  ///
  /// Same regime as `test_extreme_short_branch_dense` but exercises the sparse code path.
  /// With identical sequences, all positions are fixed (no variable positions), testing
  /// the fixed-character distribution path under near-identity transition matrices.
  #[test]
  fn test_extreme_short_branch_sparse() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let newick = "(A:1e-10,B:1e-10)root;";
    let aln = ">A\nACGT\n>B\nACGT\n";

    let (log_lh, partitions) = run_sparse_marginal_with_partitions(newick, aln, gtr)?;

    assert!(log_lh.is_finite(), "Log-likelihood is not finite: {log_lh}");
    assert!(log_lh <= 0.0, "Log-likelihood should be non-positive: {log_lh}");

    let partition = partitions[0].read_arc();
    for node_data in partition.nodes.values() {
      assert_sparse_profile_stable(&node_data.profile, 8);
    }

    Ok(())
  }

  /// Numerical stability with extremely long branches (t = 10.0), dense representation.
  ///
  /// When t is large, the transition matrix approaches the equilibrium distribution:
  /// P(t) -> pi * 1^T (each column converges to pi), because exp(lambda_i * mu * t) -> 0
  /// for all non-zero eigenvalues (which are negative for a valid rate matrix). The
  /// posterior at each position converges to the prior pi regardless of observed state.
  /// This tests whether P(t) computation handles large negative exponents without
  /// underflow.
  ///
  /// Uses completely different sequences (ACGT vs TGCA) to ensure non-trivial
  /// likelihood computation despite the near-equilibrium regime.
  #[test]
  fn test_extreme_long_branch_dense() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let newick = "(A:10.0,B:10.0)root;";
    let aln = ">A\nACGT\n>B\nTGCA\n";

    let (log_lh, partitions) = run_dense_marginal_with_partitions(newick, aln, gtr)?;

    assert!(log_lh.is_finite(), "Log-likelihood is not finite: {log_lh}");
    assert!(log_lh <= 0.0, "Log-likelihood should be non-positive: {log_lh}");

    let partition = partitions[0].read_arc();
    for node_data in partition.nodes.values() {
      assert_dense_profile_stable(&node_data.profile, 8);
    }

    Ok(())
  }

  /// Numerical stability with extremely long branches (t = 10.0), sparse representation.
  ///
  /// Same near-equilibrium regime as `test_extreme_long_branch_dense`. With completely
  /// different sequences, all positions are variable, testing the variable-position
  /// distribution path under near-equilibrium transition matrices.
  #[test]
  fn test_extreme_long_branch_sparse() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let newick = "(A:10.0,B:10.0)root;";
    let aln = ">A\nACGT\n>B\nTGCA\n";

    let (log_lh, partitions) = run_sparse_marginal_with_partitions(newick, aln, gtr)?;

    assert!(log_lh.is_finite(), "Log-likelihood is not finite: {log_lh}");
    assert!(log_lh <= 0.0, "Log-likelihood should be non-positive: {log_lh}");

    let partition = partitions[0].read_arc();
    for node_data in partition.nodes.values() {
      assert_sparse_profile_stable(&node_data.profile, 8);
    }

    Ok(())
  }

  /// Numerical stability with asymmetric branch lengths (t1 = 1e-10, t2 = 10.0), dense.
  ///
  /// One branch has a near-identity transition matrix (P(t) ~ I) while the other
  /// has a near-equilibrium matrix (P(t) -> pi * 1^T, each column converges to pi).
  /// The root must combine messages from two branches spanning ~11 orders of magnitude
  /// in transition probabilities. This tests whether message multiplication and
  /// normalization handle the scale mismatch without overflow, underflow, or loss of
  /// precision.
  #[test]
  fn test_extreme_asymmetric_branches_dense() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let newick = "(A:1e-10,B:10.0)root;";
    let aln = ">A\nACGT\n>B\nACGT\n";

    let (log_lh, partitions) = run_dense_marginal_with_partitions(newick, aln, gtr)?;

    assert!(log_lh.is_finite(), "Log-likelihood is not finite: {log_lh}");
    assert!(log_lh <= 0.0, "Log-likelihood should be non-positive: {log_lh}");

    let partition = partitions[0].read_arc();
    for node_data in partition.nodes.values() {
      assert_dense_profile_stable(&node_data.profile, 8);
    }

    Ok(())
  }

  /// Numerical stability with asymmetric branch lengths (t1 = 1e-10, t2 = 10.0), sparse.
  ///
  /// Same scale-mismatch regime as `test_extreme_asymmetric_branches_dense`. With
  /// identical sequences, all positions are fixed, testing the fixed-character path
  /// under asymmetric transition matrices.
  #[test]
  fn test_extreme_asymmetric_branches_sparse() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let newick = "(A:1e-10,B:10.0)root;";
    let aln = ">A\nACGT\n>B\nACGT\n";

    let (log_lh, partitions) = run_sparse_marginal_with_partitions(newick, aln, gtr)?;

    assert!(log_lh.is_finite(), "Log-likelihood is not finite: {log_lh}");
    assert!(log_lh <= 0.0, "Log-likelihood should be non-positive: {log_lh}");

    let partition = partitions[0].read_arc();
    for node_data in partition.nodes.values() {
      assert_sparse_profile_stable(&node_data.profile, 8);
    }

    Ok(())
  }

  // ============================================================================
  // T5: Near-zero equilibrium frequency tests
  // ============================================================================

  /// Numerical stability with skewed equilibrium frequencies and rare observed states,
  /// dense representation.
  ///
  /// Equilibrium: pi = [0.97, 0.01, 0.01, 0.01] (A dominant). Sequences use only
  /// rare states C and G (pi = 0.01 each). Small pi values stress log-space arithmetic:
  /// the site likelihood P(data_site) involves products of small transition
  /// probabilities, and the log-likelihood should be strongly negative (< -10).
  ///
  /// Near-zero pi also stresses the D^{-1/2} back-transform in the eigendecomposition
  /// (dividing by sqrt(pi_i) amplifies rounding errors when pi_i is small), testing
  /// robustness of the matrix exponential computation. The eigenvectors of the
  /// symmetrized matrix S remain orthogonal (computed via `eigh`), but the
  /// back-transformed v and v_inv matrices inherit the poor conditioning.
  #[test]
  fn test_near_zero_pi_dense() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc, true)?;
    let gtr = GTR::new(GTRParams {
      alphabet,
      mu: 1.0,
      W: None,
      pi: array![0.97, 0.01, 0.01, 0.01],
    })?;

    let newick = "(A:0.1,B:0.2)root;";
    // Use rare states to stress the small probabilities
    let aln = ">A\nCCCC\n>B\nGGGG\n";

    let (log_lh, partitions) = run_dense_marginal_with_partitions(newick, aln, gtr)?;

    assert!(log_lh.is_finite(), "Log-likelihood is not finite: {log_lh}");
    assert!(log_lh <= 0.0, "Log-likelihood should be non-positive: {log_lh}");
    // With rare states and skewed pi, log-likelihood should be very negative
    assert!(
      log_lh < -10.0,
      "Log-likelihood should be strongly negative for rare states: {log_lh}"
    );

    let partition = partitions[0].read_arc();
    for node_data in partition.nodes.values() {
      assert_dense_profile_stable(&node_data.profile, 8);
    }

    Ok(())
  }

  /// Numerical stability with skewed equilibrium frequencies and rare observed states,
  /// sparse representation.
  ///
  /// Same regime as `test_near_zero_pi_dense`. With completely different sequences
  /// (CCCC vs GGGG), all positions are variable, testing the variable-position path
  /// under skewed pi with small transition probabilities.
  #[test]
  fn test_near_zero_pi_sparse() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc, true)?;
    let gtr = GTR::new(GTRParams {
      alphabet,
      mu: 1.0,
      W: None,
      pi: array![0.97, 0.01, 0.01, 0.01],
    })?;

    let newick = "(A:0.1,B:0.2)root;";
    let aln = ">A\nCCCC\n>B\nGGGG\n";

    let (log_lh, partitions) = run_sparse_marginal_with_partitions(newick, aln, gtr)?;

    assert!(log_lh.is_finite(), "Log-likelihood is not finite: {log_lh}");
    assert!(log_lh <= 0.0, "Log-likelihood should be non-positive: {log_lh}");
    assert!(
      log_lh < -10.0,
      "Log-likelihood should be strongly negative for rare states: {log_lh}"
    );

    let partition = partitions[0].read_arc();
    for node_data in partition.nodes.values() {
      assert_sparse_profile_stable(&node_data.profile, 8);
    }

    Ok(())
  }

  /// Numerical stability with skewed pi when sequences use the dominant state, dense.
  ///
  /// Equilibrium: pi = [0.97, 0.01, 0.01, 0.01]. Sequences are all-A (the dominant
  /// state), so transition probabilities are large and the log-likelihood is much less
  /// negative than the rare-state variant. This confirms the implementation handles
  /// both tails of the probability spectrum under the same skewed model.
  #[test]
  fn test_near_zero_pi_with_dominant_state_dense() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc, true)?;
    let gtr = GTR::new(GTRParams {
      alphabet,
      mu: 1.0,
      W: None,
      pi: array![0.97, 0.01, 0.01, 0.01],
    })?;

    let newick = "(A:0.1,B:0.2)root;";
    let aln = ">A\nAAAA\n>B\nAAAA\n";

    let (log_lh, partitions) = run_dense_marginal_with_partitions(newick, aln, gtr)?;

    assert!(log_lh.is_finite(), "Log-likelihood is not finite: {log_lh}");
    assert!(log_lh <= 0.0, "Log-likelihood should be non-positive: {log_lh}");

    let partition = partitions[0].read_arc();
    for node_data in partition.nodes.values() {
      assert_dense_profile_stable(&node_data.profile, 8);
    }

    Ok(())
  }

  /// Numerical stability with extremely skewed equilibrium frequencies, dense.
  ///
  /// Equilibrium: pi = [0.9997, 0.0001, 0.0001, 0.0001]. The ratio pi_max/pi_min ~ 10^4
  /// produces a large condition number in the D^{-1/2} back-transform (dividing by
  /// sqrt(1e-4) = 0.01 amplifies errors 100x relative to the dominant state).
  /// Sequences contain all four states (ACGT), forcing the algorithm to compute
  /// likelihoods for states with pi ~ 1e-4. This stresses the back-transform
  /// conditioning and log-space arithmetic more severely than pi = [0.97, 0.01, 0.01, 0.01].
  #[test]
  fn test_extremely_skewed_pi_dense() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc, true)?;
    let gtr = GTR::new(GTRParams {
      alphabet,
      mu: 1.0,
      W: None,
      pi: array![0.9997, 0.0001, 0.0001, 0.0001],
    })?;

    let newick = "(A:0.1,B:0.1)root;";
    let aln = ">A\nACGT\n>B\nACGT\n";

    let (log_lh, partitions) = run_dense_marginal_with_partitions(newick, aln, gtr)?;

    assert!(log_lh.is_finite(), "Log-likelihood is not finite: {log_lh}");
    assert!(log_lh <= 0.0, "Log-likelihood should be non-positive: {log_lh}");

    let partition = partitions[0].read_arc();
    for node_data in partition.nodes.values() {
      assert_dense_profile_stable(&node_data.profile, 8);
    }

    Ok(())
  }

  // ============================================================================
  // T6: Rapid state transition tests (high mutation rate)
  // ============================================================================

  /// Numerical stability with high mutation rate (mu = 10), dense representation.
  ///
  /// With mu = 10 and t = 0.1, the effective evolutionary distance mu*t = 1.0 places
  /// the transition matrix in an intermediate regime where off-diagonal elements are
  /// large (high probability of state change). The rate matrix Q has off-diagonal
  /// entries Q[i,j] = W[i,j] * pi[j] (row-stochastic convention, rows sum to zero),
  /// and mu scales the eigenvalues in P(t) = v * diag(exp(lambda_i * mu * t)) * v_inv.
  /// Large mu produces large exponent magnitudes, and P(t) must remain a valid
  /// column-stochastic matrix (columns sum to 1, all entries non-negative) despite
  /// the large scaling.
  ///
  /// Uses uniform pi = [0.25, 0.25, 0.25, 0.25] (JC69-like).
  #[test]
  fn test_high_mutation_rate_dense() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc, true)?;
    let gtr = GTR::new(GTRParams {
      alphabet,
      mu: 10.0,
      W: None,
      pi: array![0.25, 0.25, 0.25, 0.25],
    })?;

    let newick = "(A:0.1,B:0.1)root;";
    let aln = ">A\nACGT\n>B\nTGCA\n";

    let (log_lh, partitions) = run_dense_marginal_with_partitions(newick, aln, gtr)?;

    assert!(log_lh.is_finite(), "Log-likelihood is not finite: {log_lh}");
    assert!(log_lh <= 0.0, "Log-likelihood should be non-positive: {log_lh}");

    let partition = partitions[0].read_arc();
    for node_data in partition.nodes.values() {
      assert_dense_profile_stable(&node_data.profile, 8);
    }

    Ok(())
  }

  /// Numerical stability with high mutation rate (mu = 10), sparse representation.
  ///
  /// Same regime as `test_high_mutation_rate_dense` (mu*t = 1.0). With completely
  /// different sequences (ACGT vs TGCA), all positions are variable.
  #[test]
  fn test_high_mutation_rate_sparse() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc, true)?;
    let gtr = GTR::new(GTRParams {
      alphabet,
      mu: 10.0,
      W: None,
      pi: array![0.25, 0.25, 0.25, 0.25],
    })?;

    let newick = "(A:0.1,B:0.1)root;";
    let aln = ">A\nACGT\n>B\nTGCA\n";

    let (log_lh, partitions) = run_sparse_marginal_with_partitions(newick, aln, gtr)?;

    assert!(log_lh.is_finite(), "Log-likelihood is not finite: {log_lh}");
    assert!(log_lh <= 0.0, "Log-likelihood should be non-positive: {log_lh}");

    let partition = partitions[0].read_arc();
    for node_data in partition.nodes.values() {
      assert_sparse_profile_stable(&node_data.profile, 8);
    }

    Ok(())
  }

  /// Numerical stability with very high mutation rate (mu = 100), dense representation.
  ///
  /// With mu = 100 and t = 0.01, the effective distance mu*t = 1.0 (same regime as
  /// mu = 10, t = 0.1). The rate matrix eigenvalues are 100x larger, producing
  /// exp(lambda_i * mu * t) values that span a wider dynamic range in the eigendecomposition.
  /// Despite the extreme scaling of Q, the product mu*t determines the physical regime,
  /// so the transition matrix should still be well-conditioned.
  ///
  /// Uses maximally different sequences (AAAA vs TTTT) with uniform pi.
  #[test]
  fn test_very_high_mutation_rate_dense() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc, true)?;
    let gtr = GTR::new(GTRParams {
      alphabet,
      mu: 100.0,
      W: None,
      pi: array![0.25, 0.25, 0.25, 0.25],
    })?;

    let newick = "(A:0.01,B:0.01)root;";
    let aln = ">A\nAAAA\n>B\nTTTT\n";

    let (log_lh, partitions) = run_dense_marginal_with_partitions(newick, aln, gtr)?;

    assert!(log_lh.is_finite(), "Log-likelihood is not finite: {log_lh}");
    assert!(log_lh <= 0.0, "Log-likelihood should be non-positive: {log_lh}");

    let partition = partitions[0].read_arc();
    for node_data in partition.nodes.values() {
      assert_dense_profile_stable(&node_data.profile, 8);
    }

    Ok(())
  }

  /// Numerical stability with high mutation rate and non-uniform equilibrium, dense.
  ///
  /// Combines mu = 10 with asymmetric pi = [0.4, 0.1, 0.2, 0.3] on a 3-taxon tree.
  /// Non-uniform pi breaks the symmetry of the JC69 model. The eigenvectors of the
  /// symmetrized matrix S remain orthogonal (computed via `eigh`), but the
  /// back-transformed v and v_inv matrices incorporate non-trivial D^{+-1/2} scaling
  /// that is absent under uniform pi. The combination of high mu and asymmetric pi
  /// tests whether the eigendecomposition and matrix exponential remain stable when
  /// neither the rate scaling nor the equilibrium is simple.
  #[test]
  fn test_high_mutation_nonuniform_pi_dense() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc, true)?;
    let gtr = GTR::new(GTRParams {
      alphabet,
      mu: 10.0,
      W: None,
      pi: array![0.4, 0.1, 0.2, 0.3],
    })?;

    let newick = "((A:0.05,B:0.05)AB:0.1,C:0.15)root;";
    let aln = ">A\nACGT\n>B\nTGCA\n>C\nGGGG\n";

    let (log_lh, partitions) = run_dense_marginal_with_partitions(newick, aln, gtr)?;

    assert!(log_lh.is_finite(), "Log-likelihood is not finite: {log_lh}");
    assert!(log_lh <= 0.0, "Log-likelihood should be non-positive: {log_lh}");

    let partition = partitions[0].read_arc();
    for node_data in partition.nodes.values() {
      assert_dense_profile_stable(&node_data.profile, 8);
    }

    Ok(())
  }

  /// Numerical stability under combined extreme conditions, dense representation.
  ///
  /// Simultaneously applies three stress factors:
  /// - Skewed pi = [0.9, 0.03, 0.04, 0.03] (poorly conditioned D^{-1/2} back-transform).
  /// - High mutation rate mu = 50 (large eigenvalue magnitudes in exp(lambda_i * mu * t)).
  /// - Very short branches t = 1e-8 (P(t) ~ I + Qt, requiring precise cancellation).
  ///
  /// The effective distance mu*t = 5e-7 is extremely small, so the transition matrix
  /// is near-identity despite the high rate. Sequences use only rare states (C, G),
  /// producing very small site likelihoods. This is the most adversarial combination
  /// for the eigendecomposition: large eigenvalues multiplied by tiny t, with
  /// poorly conditioned v/v_inv from skewed pi (dividing by sqrt(0.03) amplifies
  /// rounding errors ~5.8x), and small likelihood values from rare observations.
  #[test]
  fn test_combined_extreme_parameters_dense() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc, true)?;
    let gtr = GTR::new(GTRParams {
      alphabet,
      mu: 50.0,
      W: None,
      pi: array![0.9, 0.03, 0.04, 0.03],
    })?;

    let newick = "(A:1e-8,B:1e-8)root;";
    let aln = ">A\nCCCC\n>B\nGGGG\n";

    let (log_lh, partitions) = run_dense_marginal_with_partitions(newick, aln, gtr)?;

    assert!(log_lh.is_finite(), "Log-likelihood is not finite: {log_lh}");
    assert!(log_lh <= 0.0, "Log-likelihood should be non-positive: {log_lh}");

    let partition = partitions[0].read_arc();
    for node_data in partition.nodes.values() {
      assert_dense_profile_stable(&node_data.profile, 8);
    }

    Ok(())
  }

  /// Numerical stability under combined extreme conditions, sparse representation.
  ///
  /// Same adversarial parameter combination as `test_combined_extreme_parameters_dense`:
  /// skewed pi, high mu, very short branches, rare observed states. With completely
  /// different sequences (CCCC vs GGGG), all positions are variable, exercising the
  /// sparse variable-position path under the most extreme parameter combination.
  #[test]
  fn test_combined_extreme_parameters_sparse() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc, true)?;
    let gtr = GTR::new(GTRParams {
      alphabet,
      mu: 50.0,
      W: None,
      pi: array![0.9, 0.03, 0.04, 0.03],
    })?;

    let newick = "(A:1e-8,B:1e-8)root;";
    let aln = ">A\nCCCC\n>B\nGGGG\n";

    let (log_lh, partitions) = run_sparse_marginal_with_partitions(newick, aln, gtr)?;

    assert!(log_lh.is_finite(), "Log-likelihood is not finite: {log_lh}");
    assert!(log_lh <= 0.0, "Log-likelihood should be non-positive: {log_lh}");

    let partition = partitions[0].read_arc();
    for node_data in partition.nodes.values() {
      assert_sparse_profile_stable(&node_data.profile, 8);
    }

    Ok(())
  }
}
