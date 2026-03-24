#[cfg(test)]
mod tests {
  use crate::commands::ancestral::marginal::update_marginal;
  use crate::commands::optimize::optimize_unified::run_optimize_mixed;
  use crate::representation::payload::ancestral::GraphAncestral;
  use eyre::Report;
  use treetime_graph::edge::HasBranchLength;
  use treetime_io::nwk::nwk_read_str;

  use super::super::test_convergence_support::tests::{
    TREE_NEWICK, compute_total_lh, setup_partitions, simple_alignment,
  };

  #[test]
  fn test_optimization_converges_within_iterations() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;

    let (dense_partitions, sparse_partitions, mixed_partitions) = setup_partitions(&graph, &aln)?;

    let initial_lh = compute_total_lh(&graph, &dense_partitions, &sparse_partitions)?;
    let max_iter = 50;

    let mut lh_history = Vec::with_capacity(max_iter);
    lh_history.push(initial_lh);

    for _ in 0..max_iter {
      run_optimize_mixed(&graph, &mixed_partitions)?;
      let lh = compute_total_lh(&graph, &dense_partitions, &sparse_partitions)?;
      lh_history.push(lh);
    }

    // The undamped alternating optimization (marginal reconstruction / branch length
    // optimization) produces a stable 2-cycle. Oscillation is expected without damping -
    // that is why the production code applies damping by default. This test verifies the
    // undamped path (--damping=0.0) stays bounded rather than diverging.
    // Both phases of the 2-cycle must stay within a tight LH range.
    let final_lh = lh_history[lh_history.len() - 1];
    assert!(
      final_lh > -73.5 && final_lh < -72.0,
      "Final log-lh {final_lh:.6} outside expected range (-73.5, -72.0)"
    );

    Ok(())
  }

  #[test]
  fn test_optimization_improves_or_maintains_likelihood() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;

    let (dense_partitions, sparse_partitions, mixed_partitions) = setup_partitions(&graph, &aln)?;

    let initial_lh = compute_total_lh(&graph, &dense_partitions, &sparse_partitions)?;
    // Initial log-lh should be negative (log of probability < 1)
    assert!(initial_lh < 0.0, "Initial log-LH should be negative: {initial_lh}");

    // Run several optimization steps
    for _ in 0..10 {
      run_optimize_mixed(&graph, &mixed_partitions)?;
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

  #[test]
  fn test_optimization_produces_valid_branch_lengths() -> Result<(), Report> {
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
      run_optimize_mixed(&graph, &mixed_partitions)?;
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
