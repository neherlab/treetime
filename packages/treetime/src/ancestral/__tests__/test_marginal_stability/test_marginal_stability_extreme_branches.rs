#[cfg(test)]
mod tests {
  use super::super::test_marginal_stability_support::tests::{
    assert_dense_profile_stable, assert_sparse_profile_stable,
  };
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::test_utils::{run_dense_marginal_with_newick, run_sparse_marginal_with_newick};
  use eyre::Report;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::short_branches(       "(A:1e-10,B:1e-10)root;", ">A\nACGT\n>B\nACGT\n")]
  #[case::long_branches(        "(A:10.0,B:10.0)root;",   ">A\nACGT\n>B\nTGCA\n")]
  #[case::asymmetric_branches(  "(A:1e-10,B:10.0)root;",  ">A\nACGT\n>B\nACGT\n")]
  #[trace]
  fn test_marginal_stability_extreme_branches_dense_profiles_normalized(
    #[case] newick: &str,
    #[case] aln: &str,
  ) -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;

    let (log_lh, partitions) = run_dense_marginal_with_newick(newick, aln, &gtr)?;
    assert!(log_lh.is_finite(), "Log-likelihood is not finite: {log_lh}");
    assert!(log_lh <= 0.0, "Log-likelihood should be non-positive: {log_lh}");

    for node_data in partitions.node_states.values() {
      assert_dense_profile_stable(&node_data.profile, 8);
    }

    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::short_branches(       "(A:1e-10,B:1e-10)root;", ">A\nACGT\n>B\nACGT\n")]
  #[case::long_branches(        "(A:10.0,B:10.0)root;",   ">A\nACGT\n>B\nTGCA\n")]
  #[case::asymmetric_branches(  "(A:1e-10,B:10.0)root;",  ">A\nACGT\n>B\nACGT\n")]
  #[trace]
  fn test_marginal_stability_extreme_branches_sparse_profiles_normalized(
    #[case] newick: &str,
    #[case] aln: &str,
  ) -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;

    let (log_lh, partitions) = run_sparse_marginal_with_newick(newick, aln, &gtr)?;
    assert!(log_lh.is_finite(), "Log-likelihood is not finite: {log_lh}");
    assert!(log_lh <= 0.0, "Log-likelihood should be non-positive: {log_lh}");

    for node_data in partitions.node_states.values() {
      assert_sparse_profile_stable(&node_data.profile, 8);
    }

    Ok(())
  }
}
