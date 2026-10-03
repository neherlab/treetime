#[cfg(test)]
mod tests {
  use super::super::test_convergence_support::tests::{
    TREE_NEWICK, compute_total_lh, setup_reconstruction, simple_alignment,
  };
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::optimize::dispatch::run_optimize_mixed;
  use crate::optimize::gather::{gather_edge_contributions, gather_edge_indel_counts};
  use crate::optimize::params::BranchOptMethod;
  use eyre::Report;
  use rstest::rstest;
  use treetime_io::nwk::nwk_read_str;

  #[rustfmt::skip]
  #[rstest]
  #[case::newton(     BranchOptMethod::Newton)]
  #[case::newton_sqrt(BranchOptMethod::NewtonSqrt)]
  #[case::newton_log( BranchOptMethod::NewtonLog)]
  #[case::brent(      BranchOptMethod::Brent)]
  #[case::brent_sqrt( BranchOptMethod::BrentSqrt)]
  #[case::brent_log(  BranchOptMethod::BrentLog)]
  #[trace]
  fn test_optimization_converges_within_iterations(#[case] method: BranchOptMethod) -> Result<(), Report> {
    let aln = simple_alignment()?;
    let max_iter = 50;

    let lh_ref = {
      let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
      let graph_ref_names = nwk_parsed.names();
      let graph_ref = nwk_parsed.graph;
      let mut branch_lengths_ref = nwk_parsed.branch_lengths;
      let recon_ref = setup_reconstruction(&graph_ref, &graph_ref_names, &aln, &mut branch_lengths_ref)?;
      let total_length_ref = recon_ref.sequence_length();
      let contributions_ref = gather_edge_contributions(&graph_ref, &recon_ref)?;
      let indel_counts_ref = gather_edge_indel_counts(&graph_ref, &recon_ref);
      for _ in 0..max_iter {
        run_optimize_mixed(&graph_ref, total_length_ref, &contributions_ref, &indel_counts_ref, BranchOptMethod::BrentSqrt, &mut branch_lengths_ref)?;
      }
      let (_, lh) = compute_total_lh(&graph_ref, recon_ref, &branch_lengths_ref)?;
      lh
    };

    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;

    let reconstruction = setup_reconstruction(&graph, &names, &aln, &mut branch_lengths)?;

    let total_length = reconstruction.sequence_length();
    let contributions = gather_edge_contributions(&graph, &reconstruction)?;
    let indel_counts = gather_edge_indel_counts(&graph, &reconstruction);
    for _ in 0..max_iter {
      run_optimize_mixed(&graph, total_length, &contributions, &indel_counts, method, &mut branch_lengths)?;
    }
    let (_, final_lh) = compute_total_lh(&graph, reconstruction, &branch_lengths)?;

    let lh_diff = (final_lh - lh_ref).abs();
    assert!(
      lh_diff < 1e-2,
      "{method:?} final lh {final_lh:.6} differs from BrentSqrt reference {lh_ref:.6} by {lh_diff:.6} > 1e-2"
    );

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
  fn test_optimization_improves_or_maintains_likelihood(#[case] method: BranchOptMethod) -> Result<(), Report> {
    let aln = simple_alignment()?;
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;

    let reconstruction = setup_reconstruction(&graph, &names, &aln, &mut branch_lengths)?;

    let (reconstruction, initial_lh) = compute_total_lh(&graph, reconstruction, &branch_lengths)?;
    assert!(initial_lh < 0.0, "Initial log-LH should be negative: {initial_lh}");

    let total_length = reconstruction.sequence_length();
    let contributions = gather_edge_contributions(&graph, &reconstruction)?;
    let indel_counts = gather_edge_indel_counts(&graph, &reconstruction);
    for _ in 0..10 {
      run_optimize_mixed(&graph, total_length, &contributions, &indel_counts, method, &mut branch_lengths)?;
    }

    let (_, final_lh) = compute_total_lh(&graph, reconstruction, &branch_lengths)?;

    assert!(
      final_lh >= initial_lh,
      "Optimization regressed: {initial_lh:.6} -> {final_lh:.6}"
    );

    assert!(final_lh < 0.0, "Log-LH should be negative: {final_lh}");
    assert!(final_lh > -100.0, "Log-LH unreasonably low: {final_lh}");

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
  fn test_optimization_produces_valid_branch_lengths(#[case] method: BranchOptMethod) -> Result<(), Report> {
    let aln = simple_alignment()?;
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;

    let mut reconstruction = setup_reconstruction(&graph, &names, &aln, &mut branch_lengths)?;

    let initial_total: f64 = graph
      .get_edges()
      .filter_map(|e| branch_lengths[&e.key()])
      .sum();

    for _ in 0..10 {
      let total_length = reconstruction.sequence_length();
      let contributions = gather_edge_contributions(&graph, &reconstruction)?;
      let indel_counts = gather_edge_indel_counts(&graph, &reconstruction);
      run_optimize_mixed(&graph, total_length, &contributions, &indel_counts, method, &mut branch_lengths)?;
      (reconstruction, _) = reconstruction.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
    }

    let mut final_total = 0.0;
    for edge in graph.get_edges() {
      let branch_length = branch_lengths[&edge.key()];
      if let Some(bl) = branch_length {
        assert!(bl >= 0.0, "Branch length should be non-negative: {bl}");
        assert!(bl < 10.0, "Branch length unreasonably large: {bl}");
        final_total += bl;
      }
    }

    assert!(final_total > 0.0, "Total tree length should be positive");
    assert!(
      final_total < 50.0,
      "Total tree length unreasonably large: {final_total}"
    );

    assert!(
      final_total > initial_total * 0.1 && final_total < initial_total * 10.0,
      "Tree length changed too drastically: {initial_total} -> {final_total}"
    );

    Ok(())
  }
}
