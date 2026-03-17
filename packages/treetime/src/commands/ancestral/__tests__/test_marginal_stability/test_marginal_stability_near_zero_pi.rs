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
  fn test_near_zero_pi_dense() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let n_states = alphabet.n_canonical();
    let gtr = GTR::new(GTRParams {
      n_states,
      mu: 1.0,
      W: None,
      pi: array![0.97, 0.01, 0.01, 0.01],
    })?;

    let newick = "(A:0.1,B:0.2)root;";
    let aln = ">A\nCCCC\n>B\nGGGG\n";
    let (log_lh, partitions) = run_dense_marginal_with_partitions(newick, aln, gtr)?;

    assert!(log_lh.is_finite(), "Log-likelihood is not finite: {log_lh}");
    assert!(log_lh <= 0.0, "Log-likelihood should be non-positive: {log_lh}");
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
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let n_states = alphabet.n_canonical();
    let gtr = GTR::new(GTRParams {
      n_states,
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
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let n_states = alphabet.n_canonical();
    let gtr = GTR::new(GTRParams {
      n_states,
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
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let n_states = alphabet.n_canonical();
    let gtr = GTR::new(GTRParams {
      n_states,
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
}
