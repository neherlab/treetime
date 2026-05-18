#[cfg(test)]
mod tests {
  use crate::ancestral::marginal::update_marginal;
  use crate::optimize::params::BranchOptMethod;
  use crate::optimize::optimize_unified::run_optimize_mixed;
  use crate::representation::payload::ancestral::GraphAncestral;
  use eyre::Report;
  use rstest::rstest;
  use treetime_graph::edge::HasBranchLength;
  use treetime_io::nwk::nwk_read_str;

  use super::super::test_convergence_support::tests::{
    TREE_NEWICK, compute_total_lh, setup_partitions, simple_alignment,
  };

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

    // Reference: run BrentSqrt (the v0-matching default) on a fresh graph
    // and capture its final undamped log-likelihood. Tying the per-method
    // assertion to the reference ties this test to cross-method agreement
    // rather than to a hand-chosen LH range that would silently drift.
    let lh_ref = {
      let graph_ref: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
      let (dp_ref, sp_ref, mp_ref) = setup_partitions(&graph_ref, &aln)?;
      for _ in 0..max_iter {
        run_optimize_mixed(&graph_ref, &mp_ref, BranchOptMethod::BrentSqrt)?;
      }
      compute_total_lh(&graph_ref, &dp_ref, &sp_ref)?
    };

    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (dense_partitions, sparse_partitions, mixed_partitions) = setup_partitions(&graph, &aln)?;

    for _ in 0..max_iter {
      run_optimize_mixed(&graph, &mixed_partitions, method)?;
    }
    let final_lh = compute_total_lh(&graph, &dense_partitions, &sparse_partitions)?;

    // The undamped alternating optimization produces a stable 2-cycle. Each
    // method must converge to within 1e-2 of the BrentSqrt reference's final
    // log-likelihood. The 1e-2 tolerance accommodates the 2-cycle
    // (sub-iteration variation) for all six methods on this toy alignment.
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
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;

    let (dense_partitions, sparse_partitions, mixed_partitions) = setup_partitions(&graph, &aln)?;

    let initial_lh = compute_total_lh(&graph, &dense_partitions, &sparse_partitions)?;
    // Initial log-lh should be negative (log of probability < 1)
    assert!(initial_lh < 0.0, "Initial log-LH should be negative: {initial_lh}");

    // Run several optimization steps
    for _ in 0..10 {
      run_optimize_mixed(&graph, &mixed_partitions, method)?;
    }

    let final_lh = compute_total_lh(&graph, &dense_partitions, &sparse_partitions)?;

    // Strict non-regression: optimization must not degrade likelihood.
    // This test runs pure branch length optimization without marginal reconstruction
    // alternation, so there is no 2-cycle and likelihood should improve monotonically.
    assert!(
      final_lh >= initial_lh,
      "Optimization regressed: {initial_lh:.6} -> {final_lh:.6}"
    );

    // Final log-lh should be in reasonable range for JC69 on this alignment
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
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;

    let (dense_partitions, sparse_partitions, mixed_partitions) = setup_partitions(&graph, &aln)?;

    // Collect initial branch lengths
    let initial_total: f64 = graph
      .get_edges()
      .iter()
      .filter_map(|e| e.read_arc().payload().read_arc().branch_length())
      .sum();

    // Run several optimization iterations
    for _ in 0..10 {
      run_optimize_mixed(&graph, &mixed_partitions, method)?;
      update_marginal(&graph, &dense_partitions)?;
      update_marginal(&graph, &sparse_partitions)?;
    }

    // Verify all branch lengths are in valid range
    let mut final_total = 0.0;
    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let branch_length = edge.payload().read_arc().branch_length();
      if let Some(bl) = branch_length {
        // Branch lengths should be non-negative and reasonable (< 10 subs/site)
        assert!(bl >= 0.0, "Branch length should be non-negative: {bl}");
        assert!(bl < 10.0, "Branch length unreasonably large: {bl}");
        final_total += bl;
      }
    }

    // Total tree length should be reasonable (not zero, not huge)
    assert!(final_total > 0.0, "Total tree length should be positive");
    assert!(
      final_total < 50.0,
      "Total tree length unreasonably large: {final_total}"
    );

    // Tree length should be in same order of magnitude as initial
    // (optimization shouldn't drastically change overall scale)
    assert!(
      final_total > initial_total * 0.1 && final_total < initial_total * 10.0,
      "Tree length changed too drastically: {initial_total} -> {final_total}"
    );

    Ok(())
  }
}
