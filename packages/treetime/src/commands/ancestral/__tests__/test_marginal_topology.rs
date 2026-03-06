#[cfg(test)]
mod tests {
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::pretty_assert_ulps_eq;
  use crate::test_utils::{run_dense_marginal_with_newick, run_sparse_marginal_with_newick};
  use eyre::Report;

  // ============================================================================
  // T9: Caterpillar tree tests
  // ============================================================================

  /// Verify that dense marginal reconstruction completes on a caterpillar (maximally
  /// unbalanced) tree topology: ((((A,B),C),D),E).
  ///
  /// A caterpillar tree with N leaves has N-2 internal nodes forming a sequential spine
  /// (here 3: AB, ABC, ABCD). Each spine node except the deepest has one leaf child and
  /// one internal child; the deepest (AB) has two leaf children. Felsenstein's pruning
  /// algorithm (post-order pass) must propagate partial likelihoods through this chain,
  /// testing the sequential dependency in the backward (tips-to-root) pass.
  ///
  /// Uses identical leaf sequences (all ACGT) so all sites are invariant, testing the
  /// base case. Asserts finite, non-positive log-likelihood.
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

  /// Verify that sparse marginal reconstruction completes on a caterpillar tree.
  ///
  /// Same topology as `test_caterpillar_tree_dense_completes()` but using the sparse
  /// partition representation. With identical leaf sequences, all positions are invariant
  /// and stored as fixed distributions grouped by character. Tests that the sparse
  /// representation handles deep topologies correctly.
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

  /// Verify that dense and sparse representations produce identical log-likelihoods on a
  /// caterpillar tree with divergent sequences.
  ///
  /// The dense representation stores full probability vectors at every position, while the
  /// sparse representation stores only variable positions explicitly and groups invariant
  /// positions by character. For clean sequences (no gaps or ambiguity), both must compute
  /// the same Felsenstein likelihood. This cross-validates the two implementations on the
  /// caterpillar topology.
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

  /// Verify that dense marginal reconstruction handles a caterpillar tree where each leaf
  /// has a distinct homopolymer sequence (AAAA, CCCC, GGGG, TTTT, ACGT).
  ///
  /// Maximally divergent leaf sequences force Felsenstein's pruning to compute non-trivial
  /// partial likelihoods at every internal node, exercising the transition probability
  /// matrix multiplication P(s'|s,t) at each edge in the deep chain.
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

  /// Verify dense marginal reconstruction on a caterpillar tree with varying branch lengths
  /// spanning two orders of magnitude (0.01 to 0.8).
  ///
  /// Asymmetric branch lengths produce transition probability matrices P(t) = exp(Q*t)
  /// with very different mixing levels: short branches yield near-identity matrices while
  /// long branches approach the equilibrium distribution. This tests that the matrix
  /// exponentiation and likelihood computation remain numerically stable across this range.
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

  /// Stress test: dense marginal reconstruction on a caterpillar tree with 10 leaves (A-J).
  ///
  /// The spine has 8 internal nodes processed sequentially during the backward (post-order)
  /// pass. Tests numerical stability and stack safety when partial likelihoods are
  /// multiplied through 8 sequential Felsenstein pruning steps. With identical leaf
  /// sequences, the partial likelihood vectors should remain well-conditioned because
  /// each JC69 transition matrix P(t) = exp(Q*t) has eigenvalues 1 and exp(-4t/3),
  /// all in (0, 1] for t > 0.
  #[test]
  fn test_deep_caterpillar_tree_10_leaves() -> Result<(), Report> {
    // Create a caterpillar tree with 10 leaves (A-J)
    // Tests numerical stability and stack safety for deep recursion
    let gtr = jc69(JC69Params::default())?;
    let t = 0.1;
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

  /// Stress test: sparse marginal reconstruction on a 10-leaf caterpillar tree.
  ///
  /// Same deep topology as `test_deep_caterpillar_tree_10_leaves()` but using the sparse
  /// partition. With identical leaf sequences, all positions are invariant and handled via
  /// the fixed-distribution code path. Tests that invariant-site optimization does not
  /// introduce errors through 8 sequential spine nodes.
  #[test]
  fn test_deep_caterpillar_tree_sparse() -> Result<(), Report> {
    // Same deep caterpillar tree with sparse partition
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

  /// Verify that dense and sparse representations produce identical log-likelihoods on a
  /// 10-leaf caterpillar tree with divergent sequences.
  ///
  /// Each leaf has a distinct sequence, creating variable positions that exercise the
  /// variable-position code path in the sparse representation and the full matrix code
  /// path in the dense representation through 8 sequential spine nodes. Agreement between
  /// the two independent implementations cross-validates correctness on deep topologies.
  #[test]
  fn test_deep_tree_normalization_preserved() -> Result<(), Report> {
    // Cross-validate dense and sparse implementations on a deep caterpillar
    // tree with divergent sequences
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

  /// Stress test for numerical stability: deep caterpillar tree with branch lengths spanning
  /// over 8 orders of magnitude (1e-8 to 5.0).
  ///
  /// Very short branches (1e-8) produce JC69 transition matrices P(t) near the identity
  /// (diagonal entries ~ 1 - t), while very long branches (5.0) produce matrices near the
  /// equilibrium distribution pi = 1/4 (off-diagonal entries ~ 1/4 - 1/4*exp(-20/3)).
  /// The product of such diverse partial likelihoods through 8 sequential spine nodes
  /// tests underflow/overflow handling and the accuracy of matrix exponentiation
  /// P(t) = exp(Q*t) at extreme time scales.
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

  /// Verify Felsenstein's pruning on a star tree (polytomy) with 4 children: (A,B,C,D)root.
  ///
  /// A polytomy is a node with more than 2 children. Felsenstein's pruning generalizes
  /// from binary to arbitrary fan-out by taking the product over all children:
  ///   L_v(s) = prod_{c in children(v)} sum_{s'} P(s'|s,t_c) * L_c(s')
  /// This test verifies the product extends correctly to 4 children. All leaves share the
  /// same sequence, so the likelihood is determined purely by the tree shape and branch
  /// lengths.
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

  /// Verify dense and sparse representations produce identical log-likelihoods on a
  /// 4-child star tree (polytomy) with maximally divergent sequences.
  ///
  /// Each leaf has a different homopolymer (AAAA, CCCC, GGGG, TTTT), making every
  /// position variable. Cross-validates the two implementations on a topology where the
  /// root must integrate over 4 children simultaneously.
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

  /// Verify marginal reconstruction on a tree mixing polytomy and binary nodes:
  /// ((A,B,C)ABC,D,E)root.
  ///
  /// Internal node ABC is a 3-child polytomy, while root is also a 3-child polytomy
  /// containing one internal subtree and two leaves. Tests that Felsenstein's pruning
  /// correctly handles heterogeneous node degrees within a single tree. Dense and sparse
  /// log-likelihoods must agree.
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

  /// Verify Felsenstein's pruning on a large polytomy: a single root node with 8 leaf
  /// children (A through H).
  ///
  /// The partial likelihood at the root is the product of 8 independent edge-message
  /// terms: L_root(s) = prod_{c=A..H} sum_{s'} P(s'|s,t_c) * L_c(s'). Each leaf has a
  /// distinct 4bp sequence, but with 8 taxa and only 4 states per position, most positions
  /// have repeated states across leaves. Tests that the generalized pruning product remains
  /// numerically stable with high fan-out.
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
