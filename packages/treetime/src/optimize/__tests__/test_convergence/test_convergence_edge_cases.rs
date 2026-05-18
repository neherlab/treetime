#[cfg(test)]
mod tests {
  use crate::optimize::params::BranchOptMethod;
  use crate::optimize::dispatch::run_optimize_mixed;
  use crate::representation::payload::ancestral::GraphAncestral;
  use eyre::Report;
  use rstest::rstest;
  use treetime_graph::edge::HasBranchLength;
  use treetime_io::nwk::nwk_read_str;

  use super::super::test_convergence_support::tests::{compute_total_lh, setup_partitions, simple_alignment};

  #[rustfmt::skip]
  #[rstest]
  #[case::newton(     BranchOptMethod::Newton)]
  #[case::newton_sqrt(BranchOptMethod::NewtonSqrt)]
  #[case::newton_log( BranchOptMethod::NewtonLog)]
  #[case::brent(      BranchOptMethod::Brent)]
  #[case::brent_sqrt( BranchOptMethod::BrentSqrt)]
  #[case::brent_log(  BranchOptMethod::BrentLog)]
  #[trace]
  fn test_optimization_handles_zero_branch_lengths(#[case] method: BranchOptMethod) -> Result<(), Report> {
    // Tree with zero branch length on edge to A
    let tree_newick = "((A:0.0,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";
    let aln = simple_alignment()?;
    let graph: GraphAncestral = nwk_read_str(tree_newick)?;

    let (dense_partitions, sparse_partitions, mixed_partitions) = setup_partitions(&graph, &aln)?;

    // Run multiple optimization iterations
    for _ in 0..10 {
      run_optimize_mixed(&graph, &mixed_partitions, method)?;
    }

    let final_lh = compute_total_lh(&graph, &dense_partitions, &sparse_partitions)?;

    // Final log-lh should be in reasonable range
    assert!(final_lh < 0.0, "Log-LH should be negative: {final_lh}");
    assert!(final_lh > -100.0, "Log-LH should be reasonable: {final_lh}");

    // Branch lengths should be non-negative and bounded
    for edge in graph.get_edges() {
      let bl = edge.read_arc().payload().read_arc().branch_length();
      if let Some(bl) = bl {
        assert!(bl >= 0.0, "Branch length should be non-negative: {bl}");
        assert!(bl < 10.0, "Branch length unreasonably large: {bl}");
      }
    }

    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::newton(     BranchOptMethod::Newton)]
  #[case::newton_sqrt(BranchOptMethod::NewtonSqrt)]
  #[case::newton_log( BranchOptMethod::NewtonLog)]
  #[case::brent(      BranchOptMethod::Brent)]
  #[case::brent_sqrt( BranchOptMethod::BrentSqrt)]
  #[case::brent_log(  BranchOptMethod::BrentLog)]
  #[trace]
  fn test_optimization_handles_very_short_branches(#[case] method: BranchOptMethod) -> Result<(), Report> {
    // Tree with very short branch lengths (all 0.0001)
    let tree_newick = "((A:0.0001,B:0.0001)AB:0.0001,(C:0.0001,D:0.0001)CD:0.0001)root:0.0001;";
    let aln = simple_alignment()?;
    let graph: GraphAncestral = nwk_read_str(tree_newick)?;

    let (dense_partitions, sparse_partitions, mixed_partitions) = setup_partitions(&graph, &aln)?;

    // Run optimization iterations
    for _ in 0..10 {
      run_optimize_mixed(&graph, &mixed_partitions, method)?;
    }

    let final_lh = compute_total_lh(&graph, &dense_partitions, &sparse_partitions)?;

    // Final log-lh should be negative and reasonable
    assert!(final_lh < 0.0, "Log-LH should be negative: {final_lh}");
    assert!(final_lh > -100.0, "Log-LH should be reasonable: {final_lh}");

    // Branch lengths should be non-negative and bounded
    for edge in graph.get_edges() {
      let bl = edge.read_arc().payload().read_arc().branch_length();
      if let Some(bl) = bl {
        assert!(bl >= 0.0, "Branch length should be non-negative: {bl}");
        assert!(bl < 10.0, "Branch length unreasonably large: {bl}");
      }
    }

    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::newton(     BranchOptMethod::Newton)]
  #[case::newton_sqrt(BranchOptMethod::NewtonSqrt)]
  #[case::newton_log( BranchOptMethod::NewtonLog)]
  #[case::brent(      BranchOptMethod::Brent)]
  #[case::brent_sqrt( BranchOptMethod::BrentSqrt)]
  #[case::brent_log(  BranchOptMethod::BrentLog)]
  #[trace]
  fn test_optimization_handles_long_branches(#[case] method: BranchOptMethod) -> Result<(), Report> {
    // Tree with longer branch lengths (some > 1 sub/site)
    let tree_newick = "((A:1.0,B:2.0)AB:1.0,(C:2.0,D:1.2)CD:0.5)root:0.1;";
    let aln = simple_alignment()?;
    let graph: GraphAncestral = nwk_read_str(tree_newick)?;

    let (dense_partitions, sparse_partitions, mixed_partitions) = setup_partitions(&graph, &aln)?;

    // Run optimization iterations
    for _ in 0..10 {
      run_optimize_mixed(&graph, &mixed_partitions, method)?;
    }

    let final_lh = compute_total_lh(&graph, &dense_partitions, &sparse_partitions)?;

    // Final log-lh should be negative and reasonable
    assert!(final_lh < 0.0, "Log-LH should be negative: {final_lh}");
    assert!(final_lh > -200.0, "Log-LH should be reasonable: {final_lh}");

    // Branch lengths should be non-negative and bounded
    for edge in graph.get_edges() {
      let bl = edge.read_arc().payload().read_arc().branch_length();
      if let Some(bl) = bl {
        assert!(bl >= 0.0, "Branch length should be non-negative: {bl}");
        assert!(bl < 20.0, "Branch length unreasonably large: {bl}");
      }
    }

    Ok(())
  }
}
