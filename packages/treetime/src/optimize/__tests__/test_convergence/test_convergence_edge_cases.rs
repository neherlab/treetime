#[cfg(test)]
mod tests {
  use crate::optimize::dispatch::run_optimize_mixed;
  use crate::optimize::gather::{gather_edge_contributions, gather_edge_indel_counts, total_sequence_length};
  use crate::optimize::params::BranchOptMethod;
  use eyre::Report;
  use rstest::rstest;
  use treetime_graph::graph::Graph;
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
    let tree_newick = "((A:0.0,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";
    let aln = simple_alignment()?;
    let nwk_parsed = nwk_read_str(tree_newick)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;

    let (dense_partitions, sparse_partitions) = setup_partitions(&graph, &names, &aln, &mut branch_lengths)?;
    let total_length = total_sequence_length(&dense_partitions, &sparse_partitions);
    let contributions = gather_edge_contributions(&graph, &dense_partitions, &sparse_partitions)?;
    let indel_counts = gather_edge_indel_counts(&graph, &dense_partitions, &sparse_partitions);

    for _ in 0..10 {
      run_optimize_mixed(&graph, total_length, &contributions, &indel_counts, method, &mut branch_lengths)?;
    }

    let (dense_partitions, sparse_partitions, final_lh) = compute_total_lh(&graph, dense_partitions, sparse_partitions, &branch_lengths)?;

    assert!(final_lh < 0.0, "Log-LH should be negative: {final_lh}");
    assert!(final_lh > -100.0, "Log-LH should be reasonable: {final_lh}");

    for edge in graph.get_edges() {
      let bl = branch_lengths[&edge.key()];
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
    let tree_newick = "((A:0.0001,B:0.0001)AB:0.0001,(C:0.0001,D:0.0001)CD:0.0001)root:0.0001;";
    let aln = simple_alignment()?;
    let nwk_parsed = nwk_read_str(tree_newick)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;

    let (dense_partitions, sparse_partitions) = setup_partitions(&graph, &names, &aln, &mut branch_lengths)?;
    let total_length = total_sequence_length(&dense_partitions, &sparse_partitions);
    let contributions = gather_edge_contributions(&graph, &dense_partitions, &sparse_partitions)?;
    let indel_counts = gather_edge_indel_counts(&graph, &dense_partitions, &sparse_partitions);

    for _ in 0..10 {
      run_optimize_mixed(&graph, total_length, &contributions, &indel_counts, method, &mut branch_lengths)?;
    }

    let (dense_partitions, sparse_partitions, final_lh) = compute_total_lh(&graph, dense_partitions, sparse_partitions, &branch_lengths)?;

    assert!(final_lh < 0.0, "Log-LH should be negative: {final_lh}");
    assert!(final_lh > -100.0, "Log-LH should be reasonable: {final_lh}");

    for edge in graph.get_edges() {
      let bl = branch_lengths[&edge.key()];
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
    let tree_newick = "((A:1.0,B:2.0)AB:1.0,(C:2.0,D:1.2)CD:0.5)root:0.1;";
    let aln = simple_alignment()?;
    let nwk_parsed = nwk_read_str(tree_newick)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;

    let (dense_partitions, sparse_partitions) = setup_partitions(&graph, &names, &aln, &mut branch_lengths)?;
    let total_length = total_sequence_length(&dense_partitions, &sparse_partitions);
    let contributions = gather_edge_contributions(&graph, &dense_partitions, &sparse_partitions)?;
    let indel_counts = gather_edge_indel_counts(&graph, &dense_partitions, &sparse_partitions);

    for _ in 0..10 {
      run_optimize_mixed(&graph, total_length, &contributions, &indel_counts, method, &mut branch_lengths)?;
    }

    let (dense_partitions, sparse_partitions, final_lh) = compute_total_lh(&graph, dense_partitions, sparse_partitions, &branch_lengths)?;

    assert!(final_lh < 0.0, "Log-LH should be negative: {final_lh}");
    assert!(final_lh > -200.0, "Log-LH should be reasonable: {final_lh}");

    for edge in graph.get_edges() {
      let bl = branch_lengths[&edge.key()];
      if let Some(bl) = bl {
        assert!(bl >= 0.0, "Branch length should be non-negative: {bl}");
        assert!(bl < 20.0, "Branch length unreasonably large: {bl}");
      }
    }

    Ok(())
  }
}
