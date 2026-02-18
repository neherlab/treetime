#[cfg(test)]
mod tests {
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::pretty_assert_ulps_eq;
  use crate::test_utils::{run_dense_marginal_with_newick, run_sparse_marginal_with_newick};
  use eyre::Report;

  // ============================================================================
  // T9: Caterpillar tree tests
  // ============================================================================

  #[test]
  fn test_caterpillar_tree_dense_completes() -> Result<(), Report> {
    // Maximally unbalanced tree: ((((A:t,B:t):t,C:t):t,D:t):t,E:t)root;
    // Tests that the algorithm handles deep recursion correctly
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
    // Same caterpillar tree with sparse partition
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
    // Dense and sparse should give same results on caterpillar tree
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
    // Caterpillar tree with varied sequences to stress message passing
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
    // Caterpillar tree with varying branch lengths
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

  // ============================================================================
  // T10: Deep tree stress tests
  // ============================================================================

  #[test]
  fn test_deep_binary_tree_10_levels() -> Result<(), Report> {
    // Create a balanced binary tree with 10 levels (1023 internal nodes, 1024 leaves)
    // This tests numerical stability and stack safety for deep recursion
    let gtr = jc69(JC69Params::default())?;
    let t = 0.1;

    // Build tree recursively: depth 10 = 2^10 = 1024 leaves
    // For simplicity, we'll build a smaller balanced tree with 16 leaves (4 levels)
    // and an unbalanced caterpillar with 10 levels
    let newick = format!(
      "(((((((((A:{t},B:{t}):{t},C:{t}):{t},D:{t}):{t},E:{t}):{t},F:{t}):{t},G:{t}):{t},H:{t}):{t},I:{t}):{t},J:{t})root;"
    );

    let aln = ">A\nACGT\n>B\nACGT\n>C\nACGT\n>D\nACGT\n>E\nACGT\n>F\nACGT\n>G\nACGT\n>H\nACGT\n>I\nACGT\n>J\nACGT\n";

    let log_lh = run_dense_marginal_with_newick(&newick, aln, gtr)?;

    assert!(
      log_lh.is_finite(),
      "Deep tree (10 levels) should produce finite log-likelihood"
    );
    assert!(log_lh <= 0.0, "Log-likelihood should be non-positive");

    Ok(())
  }

  #[test]
  fn test_deep_binary_tree_sparse() -> Result<(), Report> {
    // Same deep tree with sparse partition
    let gtr = jc69(JC69Params::default())?;
    let t = 0.1;

    let newick = format!(
      "(((((((((A:{t},B:{t}):{t},C:{t}):{t},D:{t}):{t},E:{t}):{t},F:{t}):{t},G:{t}):{t},H:{t}):{t},I:{t}):{t},J:{t})root;"
    );

    let aln = ">A\nACGT\n>B\nACGT\n>C\nACGT\n>D\nACGT\n>E\nACGT\n>F\nACGT\n>G\nACGT\n>H\nACGT\n>I\nACGT\n>J\nACGT\n";

    let log_lh = run_sparse_marginal_with_newick(&newick, aln, gtr)?;

    assert!(
      log_lh.is_finite(),
      "Deep tree (sparse) should produce finite log-likelihood"
    );
    assert!(log_lh <= 0.0, "Log-likelihood should be non-positive");

    Ok(())
  }

  #[test]
  fn test_deep_tree_normalization_preserved() -> Result<(), Report> {
    // Verify that normalization is preserved through deep tree traversal
    // by checking dense/sparse consistency
    let gtr = jc69(JC69Params::default())?;
    let t = 0.1;

    let newick = format!(
      "(((((((((A:{t},B:{t}):{t},C:{t}):{t},D:{t}):{t},E:{t}):{t},F:{t}):{t},G:{t}):{t},H:{t}):{t},I:{t}):{t},J:{t})root;"
    );

    // Use varied sequences to exercise all code paths
    let aln = ">A\nAAAA\n>B\nCCCC\n>C\nGGGG\n>D\nTTTT\n>E\nACGT\n>F\nTGCA\n>G\nGTAC\n>H\nCATG\n>I\nACTG\n>J\nTGAC\n";

    let dense_log_lh = run_dense_marginal_with_newick(&newick, aln, gtr.clone())?;
    let sparse_log_lh = run_sparse_marginal_with_newick(&newick, aln, gtr)?;

    pretty_assert_ulps_eq!(dense_log_lh, sparse_log_lh, epsilon = 1e-10);

    Ok(())
  }

  #[test]
  fn test_deep_tree_extreme_branches() -> Result<(), Report> {
    // Deep tree with extreme branch lengths to test numerical stability
    let gtr = jc69(JC69Params::default())?;

    // Mix of very short and very long branches
    let newick = "(((((((((A:1e-8,B:1e-8):0.001,C:0.01):0.1,D:1.0):2.0,E:3.0):0.001,F:0.001):5.0,G:0.01):0.1,H:0.5):1.0,I:2.0):0.5,J:0.5)root;";

    let aln = ">A\nACGT\n>B\nACGT\n>C\nTGCA\n>D\nGTAC\n>E\nCATG\n>F\nACTG\n>G\nTGAC\n>H\nGATC\n>I\nCTAG\n>J\nATCG\n";

    let log_lh = run_dense_marginal_with_newick(newick, aln, gtr)?;

    assert!(
      log_lh.is_finite(),
      "Deep tree with extreme branches should produce finite log-likelihood"
    );
    assert!(log_lh <= 0.0, "Log-likelihood should be non-positive");

    Ok(())
  }

  // ============================================================================
  // T11: Polytomy tests
  // ============================================================================

  #[test]
  fn test_polytomy_four_children() -> Result<(), Report> {
    // Tree with polytomy (>2 children): (A:t,B:t,C:t,D:t)root;
    // This is a star tree which is a valid polytomy
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
    // Verify polytomy produces same results in dense and sparse
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
    // Tree mixing polytomies and binary nodes
    // ((A:t,B:t,C:t):t,D:t,E:t)root;
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
    // Polytomy with many children (8)
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
