#[cfg(test)]
mod tests {
  use super::super::test_marginal_stability_support::tests::{
    assert_dense_profile_stable, assert_sparse_profile_stable, run_dense_marginal_with_partitions,
    run_sparse_marginal_with_partitions,
  };
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::gtr::gtr::{GTR, GTRParams};
  use eyre::Report;
  use ndarray::array;

  #[test]
  fn test_high_mutation_rate_dense() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let n_states = alphabet.n_canonical();
    let gtr = GTR::new(GTRParams {
      n_states,
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
    for node_data in partition.data.nodes.values() {
      assert_dense_profile_stable(&node_data.profile, 8);
    }

    Ok(())
  }

  #[test]
  fn test_high_mutation_rate_sparse() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let n_states = alphabet.n_canonical();
    let gtr = GTR::new(GTRParams {
      n_states,
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
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let n_states = alphabet.n_canonical();
    let gtr = GTR::new(GTRParams {
      n_states,
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
    for node_data in partition.data.nodes.values() {
      assert_dense_profile_stable(&node_data.profile, 8);
    }

    Ok(())
  }

  #[test]
  fn test_high_mutation_nonuniform_pi_dense() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let n_states = alphabet.n_canonical();
    let gtr = GTR::new(GTRParams {
      n_states,
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
    for node_data in partition.data.nodes.values() {
      assert_dense_profile_stable(&node_data.profile, 8);
    }

    Ok(())
  }

  #[test]
  fn test_combined_extreme_parameters_dense() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let n_states = alphabet.n_canonical();
    let gtr = GTR::new(GTRParams {
      n_states,
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
    for node_data in partition.data.nodes.values() {
      assert_dense_profile_stable(&node_data.profile, 8);
    }

    Ok(())
  }

  #[test]
  fn test_combined_extreme_parameters_sparse() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let n_states = alphabet.n_canonical();
    let gtr = GTR::new(GTRParams {
      n_states,
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
