#[cfg(test)]
mod tests {
  use super::super::test_marginal_stability_support::tests::{
    assert_dense_profile_stable, assert_sparse_profile_stable, run_dense_marginal_with_partitions,
    run_sparse_marginal_with_partitions,
  };
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use eyre::Report;

  #[test]
  fn test_extreme_short_branch_dense() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let newick = "(A:1e-10,B:1e-10)root;";
    let aln = ">A\nACGT\n>B\nACGT\n";

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
    let gtr = jc69(JC69Params::default())?;
    let newick = "(A:10.0,B:10.0)root;";
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
    let gtr = jc69(JC69Params::default())?;
    let newick = "(A:1e-10,B:10.0)root;";
    let aln = ">A\nACGT\n>B\nACGT\n";

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
}
