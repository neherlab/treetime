#[cfg(test)]
mod tests {
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::pretty_assert_ulps_eq;
  use crate::test_utils::{run_dense_marginal_with_newick, run_sparse_marginal_with_newick};
  use eyre::Report;
  use rstest::rstest;

  #[rstest]
  #[case::identical_sequences(
    ">A\nACGT\n>B\nACGT\n>C\nACGT\n>D\nACGT\n>E\nACGT\n>F\nACGT\n>G\nACGT\n>H\nACGT\n>I\nACGT\n>J\nACGT\n"
  )]
  #[case::diverse_sequences(
    ">A\nAAAA\n>B\nCCCC\n>C\nGGGG\n>D\nTTTT\n>E\nACGT\n>F\nTGCA\n>G\nGTAC\n>H\nCATG\n>I\nACTG\n>J\nTGAC\n"
  )]
  #[trace]
  fn test_marginal_topology_deep_caterpillar_dense_matches_sparse(#[case] aln: &str) -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let t = 0.1;
    let newick = format!(
      "(((((((((A:{t},B:{t}):{t},C:{t}):{t},D:{t}):{t},E:{t}):{t},F:{t}):{t},G:{t}):{t},H:{t}):{t},I:{t}):{t},J:{t})root;"
    );

    let (dense_log_lh, _) = run_dense_marginal_with_newick(&newick, aln, &gtr)?;
    let (sparse_log_lh, _) = run_sparse_marginal_with_newick(&newick, aln, &gtr)?;

    assert!(dense_log_lh.is_finite(), "log-likelihood is not finite: {dense_log_lh}");
    assert!(dense_log_lh <= 0.0, "log-likelihood is positive: {dense_log_lh}");
    pretty_assert_ulps_eq!(dense_log_lh, sparse_log_lh, epsilon = 1e-10);
    Ok(())
  }

  #[test]
  fn test_marginal_topology_deep_tree_extreme_branches_log_lh_finite() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let newick = "(((((((((A:1e-8,B:1e-8):0.001,C:0.01):0.1,D:1.0):2.0,E:3.0):0.001,F:0.001):5.0,G:0.01):0.1,H:0.5):1.0,I:2.0):0.5,J:0.5)root;";
    let aln = ">A\nACGT\n>B\nACGT\n>C\nTGCA\n>D\nGTAC\n>E\nCATG\n>F\nACTG\n>G\nTGAC\n>H\nGATC\n>I\nCTAG\n>J\nATCG\n";

    let (log_lh, _) = run_dense_marginal_with_newick(newick, aln, &gtr)?;
    assert!(
      log_lh.is_finite(),
      "Deep tree with extreme branches should produce finite log-likelihood"
    );
    assert!(log_lh <= 0.0, "Log-likelihood should be non-positive");

    Ok(())
  }
}
