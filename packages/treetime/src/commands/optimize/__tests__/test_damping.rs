#[cfg(test)]
mod tests {
  use crate::commands::ancestral::marginal::update_marginal;
  use crate::commands::optimize::__tests__::test_convergence::test_convergence_support::tests::{
    TREE_NEWICK, compute_total_lh, setup_partitions, simple_alignment,
  };
  use crate::commands::optimize::optimize_unified::run_optimize_mixed;
  use crate::commands::optimize::run::{apply_damping, save_branch_lengths};
  use crate::representation::payload::ancestral::GraphAncestral;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use itertools::izip;
  use rstest::rstest;
  use treetime_graph::edge::HasBranchLength;
  use treetime_io::nwk::nwk_read_str;

  #[test]
  fn test_save_branch_lengths_captures_all_edges() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let saved = save_branch_lengths(&graph);
    let edges = graph.get_edges();
    assert_eq!(saved.len(), edges.len());
    for (bl, edge_ref) in izip!(&saved, &edges) {
      let expected = edge_ref.read_arc().payload().read_arc().branch_length().unwrap_or(0.0);
      assert_abs_diff_eq!(*bl, expected, epsilon = 1e-15);
    }
    Ok(())
  }

  #[test]
  fn test_apply_damping_zero_is_noop() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let original = save_branch_lengths(&graph);

    // Manually change branch lengths to simulate optimization
    for edge_ref in graph.get_edges() {
      let mut edge = edge_ref.write_arc().payload().write_arc();
      let bl = edge.branch_length().unwrap_or(0.0);
      edge.set_branch_length(Some(bl * 2.0));
    }

    let after_optim = save_branch_lengths(&graph);
    apply_damping(&graph, &original, 0.0, 0);
    let after_damping = save_branch_lengths(&graph);

    // damping=0.0 should leave the optimized values untouched
    for (damped, optimized) in izip!(&after_damping, &after_optim) {
      assert_abs_diff_eq!(*damped, *optimized, epsilon = 1e-15);
    }
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::iter_0( 0, 0.750)]
  #[case::iter_1( 1, 0.5625)]
  #[case::iter_2( 2, 0.421875)]
  #[case::iter_4( 4, 0.2373046875)]
  #[case::iter_9( 9, 0.056313514709472656)]
  #[trace]
  fn test_apply_damping_weights_match_v0(#[case] iteration: usize, #[case] expected_old_weight: f64) -> Result<(), Report> {
    let damping = 0.75;
    let graph: GraphAncestral = nwk_read_str("(A:1.0,B:1.0)root:0.0;")?;
    let old_bls = save_branch_lengths(&graph);

    // Set all branch lengths to a known "optimized" value
    for edge_ref in graph.get_edges() {
      edge_ref.write_arc().payload().write_arc().set_branch_length(Some(0.0));
    }

    apply_damping(&graph, &old_bls, damping, iteration);

    // bl = 0.0 * (1 - old_weight) + 1.0 * old_weight = old_weight
    for edge_ref in graph.get_edges() {
      let bl = edge_ref.read_arc().payload().read_arc().branch_length().unwrap_or(0.0);
      assert_abs_diff_eq!(bl, expected_old_weight, epsilon = 1e-15);
    }
    Ok(())
  }

  #[test]
  fn test_apply_damping_blends_correctly() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("(A:0.1,B:0.2)root:0.0;")?;
    let old_bls = save_branch_lengths(&graph);

    // Set "optimized" branch lengths
    for edge_ref in graph.get_edges() {
      let mut edge = edge_ref.write_arc().payload().write_arc();
      let bl = edge.branch_length().unwrap_or(0.0);
      edge.set_branch_length(Some(bl * 3.0));
    }

    apply_damping(&graph, &old_bls, 0.75, 0);

    // At iteration 0, damping_factor = 0.75, new_weight = 0.25
    // Edge A: 0.3 * 0.25 + 0.1 * 0.75 = 0.075 + 0.075 = 0.15
    // Edge B: 0.6 * 0.25 + 0.2 * 0.75 = 0.15 + 0.15 = 0.30
    // Epsilon accounts for Newick float parsing roundtrip
    let damped = save_branch_lengths(&graph);
    assert_abs_diff_eq!(damped[0], 0.15, epsilon = 1e-8);
    assert_abs_diff_eq!(damped[1], 0.30, epsilon = 1e-8);
    Ok(())
  }

  #[test]
  fn test_apply_damping_new_weight_increases_with_iteration() -> Result<(), Report> {
    let damping = 0.75;
    let old_bl = 1.0;
    let optimized_bl = 2.0;

    let mut prev_damped = old_bl;
    for iteration in 0..10 {
      let graph: GraphAncestral = nwk_read_str("(A:1.0)root:0.0;")?;
      let old_bls = save_branch_lengths(&graph);
      graph
        .get_edges()
        .first()
        .unwrap()
        .write_arc()
        .payload()
        .write_arc()
        .set_branch_length(Some(optimized_bl));

      apply_damping(&graph, &old_bls, damping, iteration);

      let damped = graph
        .get_edges()
        .first()
        .unwrap()
        .read_arc()
        .payload()
        .read_arc()
        .branch_length()
        .unwrap_or(0.0);

      // Each subsequent iteration should give more weight to the optimized value
      assert!(
        damped > prev_damped || iteration == 0,
        "Iteration {iteration}: damped {damped} should be > previous {prev_damped}"
      );
      prev_damped = damped;
    }
    Ok(())
  }

  #[test]
  fn test_damped_optimization_converges() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (dense_partitions, sparse_partitions, mixed_partitions) = setup_partitions(&graph, &aln)?;

    let max_iter = 20;
    let damping = 0.75;
    let mut lh_history = Vec::with_capacity(max_iter);

    for i in 0..max_iter {
      let lh = compute_total_lh(&graph, &dense_partitions, &sparse_partitions)?;
      lh_history.push(lh);

      let old_bls = save_branch_lengths(&graph);
      run_optimize_mixed(&graph, &mixed_partitions)?;
      apply_damping(&graph, &old_bls, damping, i);
    }

    let lh = compute_total_lh(&graph, &dense_partitions, &sparse_partitions)?;
    lh_history.push(lh);

    // Variance over last 5 iterations should be small (convergence)
    let last_5: Vec<f64> = lh_history.iter().rev().take(5).copied().collect();
    let mean = last_5.iter().sum::<f64>() / 5.0;
    let variance = last_5.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / 5.0;
    assert!(
      variance < 1.0,
      "Damped optimization should converge: variance of last 5 iterations = {variance}"
    );

    // Final likelihood should be in reasonable range
    let final_lh = *lh_history.last().unwrap();
    assert!(
      final_lh > -100.0 && final_lh < -30.0,
      "Final log-lh {final_lh} outside expected range [-100, -30]"
    );

    Ok(())
  }

  #[test]
  fn test_damped_optimization_does_not_regress_significantly() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (dense_partitions, sparse_partitions, mixed_partitions) = setup_partitions(&graph, &aln)?;

    let initial_lh = compute_total_lh(&graph, &dense_partitions, &sparse_partitions)?;

    let damping = 0.75;
    for i in 0..10 {
      let old_bls = save_branch_lengths(&graph);
      run_optimize_mixed(&graph, &mixed_partitions)?;
      apply_damping(&graph, &old_bls, damping, i);
      update_marginal(&graph, &dense_partitions)?;
      update_marginal(&graph, &sparse_partitions)?;
    }

    let final_lh = compute_total_lh(&graph, &dense_partitions, &sparse_partitions)?;
    assert!(
      final_lh >= initial_lh - 1.0,
      "Damped optimization should not regress significantly: {initial_lh} -> {final_lh}"
    );
    Ok(())
  }
}
