#[cfg(test)]
mod tests {
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::pretty_assert_ulps_eq;
  use crate::test_utils::{run_dense_marginal_with_newick, run_sparse_marginal_with_newick};
  use eyre::Report;

  #[test]
  fn test_caterpillar_tree_dense_completes() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let t = 0.1;

    let newick = format!("((((A:{t},B:{t})AB:{t},C:{t})ABC:{t},D:{t})ABCD:{t},E:{t})root;");
    let aln = ">A\nACGT\n>B\nACGT\n>C\nACGT\n>D\nACGT\n>E\nACGT\n";

    let log_lh = run_dense_marginal_with_newick(&newick, aln, gtr)?;
    assert!(
      log_lh.is_finite(),
      "Caterpillar tree should produce finite log-likelihood"
    );
    assert!(log_lh <= 0.0, "Log-likelihood should be non-positive");

    Ok(())
  }

  #[test]
  fn test_caterpillar_tree_sparse_completes() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let t = 0.1;

    let newick = format!("((((A:{t},B:{t})AB:{t},C:{t})ABC:{t},D:{t})ABCD:{t},E:{t})root;");
    let aln = ">A\nACGT\n>B\nACGT\n>C\nACGT\n>D\nACGT\n>E\nACGT\n";

    let log_lh = run_sparse_marginal_with_newick(&newick, aln, gtr)?;
    assert!(
      log_lh.is_finite(),
      "Sparse caterpillar tree should produce finite log-likelihood"
    );
    assert!(log_lh <= 0.0, "Log-likelihood should be non-positive");

    Ok(())
  }

  #[test]
  fn test_caterpillar_tree_dense_sparse_consistency() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let t = 0.1;

    let newick = format!("((((A:{t},B:{t})AB:{t},C:{t})ABC:{t},D:{t})ABCD:{t},E:{t})root;");
    let aln = ">A\nACGTACGT\n>B\nTGCATGCA\n>C\nGTACGTAC\n>D\nCATGCATG\n>E\nACGTTGCA\n";

    let dense_log_lh = run_dense_marginal_with_newick(&newick, aln, gtr.clone())?;
    let sparse_log_lh = run_sparse_marginal_with_newick(&newick, aln, gtr)?;

    pretty_assert_ulps_eq!(dense_log_lh, sparse_log_lh, epsilon = 1e-10);
    Ok(())
  }

  #[test]
  fn test_caterpillar_tree_varied_sequences() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let t = 0.1;

    let newick = format!("((((A:{t},B:{t})AB:{t},C:{t})ABC:{t},D:{t})ABCD:{t},E:{t})root;");
    let aln = ">A\nAAAA\n>B\nCCCC\n>C\nGGGG\n>D\nTTTT\n>E\nACGT\n";

    let log_lh = run_dense_marginal_with_newick(&newick, aln, gtr)?;
    assert!(
      log_lh.is_finite(),
      "Varied sequences should produce finite log-likelihood"
    );
    assert!(log_lh <= 0.0, "Log-likelihood should be non-positive");

    Ok(())
  }

  #[test]
  fn test_caterpillar_tree_asymmetric_branches() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let newick = "((((A:0.01,B:0.02)AB:0.05,C:0.1)ABC:0.2,D:0.3)ABCD:0.5,E:0.8)root;";
    let aln = ">A\nACGT\n>B\nACGT\n>C\nTGCA\n>D\nGTAC\n>E\nCATG\n";

    let log_lh = run_dense_marginal_with_newick(newick, aln, gtr)?;
    assert!(
      log_lh.is_finite(),
      "Asymmetric branches should produce finite log-likelihood"
    );
    assert!(log_lh <= 0.0, "Log-likelihood should be non-positive");

    Ok(())
  }
}
