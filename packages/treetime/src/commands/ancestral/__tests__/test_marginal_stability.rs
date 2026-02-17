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

  #[test]
  fn test_extreme_short_branch_dense() -> Result<(), Report> {
    // Very short branch: t = 1e-10 (nearly identical sequences expected)
    // Verifies no numerical instability from near-identity transition matrices
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

  #[test]
  fn test_extreme_long_branch_dense() -> Result<(), Report> {
    // Very long branch: t = 10.0 (near equilibrium)
    // Verifies transition matrix computation doesn't overflow/underflow
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

  #[test]
  fn test_extreme_asymmetric_branches_dense() -> Result<(), Report> {
    // Mix of very short and very long branches
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

  #[test]
  fn test_near_zero_pi_dense() -> Result<(), Report> {
    // GTR with very skewed equilibrium frequencies
    // pi = [0.97, 0.01, 0.01, 0.01] - 'A' dominant
    // Verifies fixed state distributions don't underflow
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

  #[test]
  fn test_near_zero_pi_with_dominant_state_dense() -> Result<(), Report> {
    // Sequences use the dominant state - should have higher likelihood
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

  #[test]
  fn test_extremely_skewed_pi_dense() -> Result<(), Report> {
    // Even more extreme: one frequency very close to 1
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

  #[test]
  fn test_high_mutation_rate_dense() -> Result<(), Report> {
    // High mutation rate (mu = 10.0) with short branch (t = 0.1)
    // gives expQt with large off-diagonal elements
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

  #[test]
  fn test_very_high_mutation_rate_dense() -> Result<(), Report> {
    // Even higher mutation rate (mu = 100.0)
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

  #[test]
  fn test_high_mutation_nonuniform_pi_dense() -> Result<(), Report> {
    // High mutation rate with non-uniform equilibrium
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

  #[test]
  fn test_combined_extreme_parameters_dense() -> Result<(), Report> {
    // Combine multiple extreme conditions:
    // - Skewed pi
    // - High mutation rate
    // - Very short branch
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
