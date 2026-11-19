#[cfg(test)]
mod tests {
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::pretty_assert_ulps_eq;
  use crate::test_utils::{run_dense_marginal_with_newick, run_sparse_marginal_with_newick};
  use eyre::Report;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::four_children_identical(
    "(A:0.1,B:0.1,C:0.1,D:0.1)root;",
    ">A\nACGT\n>B\nACGT\n>C\nACGT\n>D\nACGT\n",
  )]
  #[case::four_children_diverse(
    "(A:0.1,B:0.1,C:0.1,D:0.1)root;",
    ">A\nAAAA\n>B\nCCCC\n>C\nGGGG\n>D\nTTTT\n",
  )]
  #[case::polytomy_with_nested_polytomy(
    "((A:0.1,B:0.1,C:0.1)ABC:0.1,D:0.1,E:0.1)root;",
    ">A\nACGT\n>B\nTGCA\n>C\nGTAC\n>D\nCATG\n>E\nACTG\n",
  )]
  #[case::eight_children(
    "(A:0.1,B:0.1,C:0.1,D:0.1,E:0.1,F:0.1,G:0.1,H:0.1)root;",
    ">A\nACGT\n>B\nTGCA\n>C\nGTAC\n>D\nCATG\n>E\nACTG\n>F\nTGAC\n>G\nGATC\n>H\nCTAG\n",
  )]
  #[trace]
  fn test_marginal_topology_polytomy_dense_matches_sparse(#[case] newick: &str, #[case] aln: &str) -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;

    let (dense_log_lh, _) = run_dense_marginal_with_newick(newick, aln, &gtr)?;
    let (sparse_log_lh, _) = run_sparse_marginal_with_newick(newick, aln, &gtr)?;

    assert!(dense_log_lh.is_finite(), "log-likelihood is not finite: {dense_log_lh}");
    assert!(dense_log_lh <= 0.0, "log-likelihood is positive: {dense_log_lh}");
    pretty_assert_ulps_eq!(dense_log_lh, sparse_log_lh, epsilon = 1e-10);
    Ok(())
  }
}
