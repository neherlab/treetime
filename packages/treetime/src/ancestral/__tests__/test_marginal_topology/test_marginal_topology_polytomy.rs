#[cfg(test)]
mod tests {
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::pretty_assert_ulps_eq;
  use crate::test_utils::{run_dense_marginal_with_newick, run_sparse_marginal_with_newick};
  use eyre::Report;

  #[test]
  fn test_polytomy_four_children() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let t = 0.1;

    let newick = format!("(A:{t},B:{t},C:{t},D:{t})root;");
    let aln = ">A\nACGT\n>B\nACGT\n>C\nACGT\n>D\nACGT\n";

    let log_lh = run_dense_marginal_with_newick(&newick, aln, gtr)?;
    assert!(log_lh.is_finite(), "Polytomy should produce finite log-likelihood");
    assert!(log_lh <= 0.0, "Log-likelihood should be non-positive");

    Ok(())
  }

  #[test]
  fn test_polytomy_dense_sparse_consistency() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let t = 0.1;

    let newick = format!("(A:{t},B:{t},C:{t},D:{t})root;");
    let aln = ">A\nAAAA\n>B\nCCCC\n>C\nGGGG\n>D\nTTTT\n";

    let dense_log_lh = run_dense_marginal_with_newick(&newick, aln, gtr.clone())?;
    let sparse_log_lh = run_sparse_marginal_with_newick(&newick, aln, gtr)?;

    pretty_assert_ulps_eq!(dense_log_lh, sparse_log_lh, epsilon = 1e-10);
    Ok(())
  }

  #[test]
  fn test_mixed_polytomy_binary() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let t = 0.1;

    let newick = format!("((A:{t},B:{t},C:{t})ABC:{t},D:{t},E:{t})root;");
    let aln = ">A\nACGT\n>B\nTGCA\n>C\nGTAC\n>D\nCATG\n>E\nACTG\n";

    let dense_log_lh = run_dense_marginal_with_newick(&newick, aln, gtr.clone())?;
    let sparse_log_lh = run_sparse_marginal_with_newick(&newick, aln, gtr)?;

    assert!(
      dense_log_lh.is_finite(),
      "Mixed tree should produce finite log-likelihood"
    );
    pretty_assert_ulps_eq!(dense_log_lh, sparse_log_lh, epsilon = 1e-10);

    Ok(())
  }

  #[test]
  fn test_large_polytomy() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let t = 0.1;

    let newick = format!("(A:{t},B:{t},C:{t},D:{t},E:{t},F:{t},G:{t},H:{t})root;");
    let aln = ">A\nACGT\n>B\nTGCA\n>C\nGTAC\n>D\nCATG\n>E\nACTG\n>F\nTGAC\n>G\nGATC\n>H\nCTAG\n";

    let log_lh = run_dense_marginal_with_newick(&newick, aln, gtr)?;
    assert!(
      log_lh.is_finite(),
      "Large polytomy should produce finite log-likelihood"
    );
    assert!(log_lh <= 0.0, "Log-likelihood should be non-positive");

    Ok(())
  }
}
