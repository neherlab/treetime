#[cfg(test)]
mod tests {
  use crate::optimize::__tests__::test_convergence::test_convergence_support::tests::{
    TREE_NEWICK, compute_total_lh, setup_partitions, simple_alignment,
  };
  use crate::optimize::iteration::apply_damping;
  use crate::optimize::params::{BranchOptMethod, TopologyOps};
  use crate::optimize::run_loop::run_optimize_loop;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_io::nwk::nwk_read_str;

  #[test]
  fn test_parse_branch_lengths_covers_all_edges() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let saved = branch_lengths;
    let edges = graph.get_edges().collect::<Vec<_>>();
    assert_eq!(saved.len(), edges.len());
    for edge_ref in &edges {
      let edge = edge_ref;
      assert!(saved.contains_key(&edge.key()));
    }
    Ok(())
  }

  #[test]
  fn test_apply_damping_zero_is_noop() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let original = branch_lengths;

    let mut optimized: BTreeMap<GraphEdgeKey, Option<f64>> =
      original.iter().map(|(&key, &bl)| (key, bl.map(|b| b * 2.0))).collect();
    let after_optim = optimized.clone();

    apply_damping(&mut optimized, &original, 0.0, 0);

    assert_eq!(optimized, after_optim);
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
    let nwk_parsed = nwk_read_str("(A:1.0,B:1.0)root:0.0;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let old_bls = branch_lengths;

    let mut bls: BTreeMap<GraphEdgeKey, Option<f64>> = old_bls.keys().map(|&key| (key, Some(0.0))).collect();

    apply_damping(&mut bls, &old_bls, damping, iteration);

    for bl in bls.values() {
      assert_abs_diff_eq!(bl.unwrap(), expected_old_weight, epsilon = 1e-15);
    }
    Ok(())
  }

  #[test]
  fn test_apply_damping_blends_correctly() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:0.1,B:0.2)root:0.0;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let old_bls = branch_lengths;

    let mut bls: BTreeMap<GraphEdgeKey, Option<f64>> =
      old_bls.iter().map(|(&key, &bl)| (key, bl.map(|b| b * 3.0))).collect();

    apply_damping(&mut bls, &old_bls, 0.75, 0);

    for (&key, &old) in &old_bls {
      let damped = bls[&key].unwrap();
      let expected = if (old.unwrap() - 0.1).abs() < 1e-9 { 0.15 } else { 0.30 };
      assert_abs_diff_eq!(damped, expected, epsilon = 1e-8);
    }
    Ok(())
  }

  #[test]
  fn test_apply_damping_new_weight_increases_with_iteration() -> Result<(), Report> {
    let damping = 0.75;
    let old_bl = 1.0;
    let optimized_bl = 2.0;

    let mut prev_damped = old_bl;
    for iteration in 0..10 {
      let nwk_parsed = nwk_read_str("(A:1.0)root:0.0;")?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;
      let graph: Graph = graph;
      let old_bls = branch_lengths;
      let mut bls: BTreeMap<GraphEdgeKey, Option<f64>> = old_bls.keys().map(|&key| (key, Some(optimized_bl))).collect();

      apply_damping(&mut bls, &old_bls, damping, iteration);

      let damped = bls.values().next().unwrap().unwrap();

      assert!(
        damped > prev_damped || iteration == 0,
        "Iteration {iteration}: damped {damped} should be > previous {prev_damped}"
      );
      prev_damped = damped;
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
  fn test_damped_optimization_converges(#[case] method: BranchOptMethod) -> Result<(), Report> {
    let aln = simple_alignment()?;
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let (dense_partitions, sparse_partitions) = setup_partitions(&graph, &names, &aln, &mut branch_lengths)?;

    let max_iter = 10;
    let damping = 0.75;
    let dp = 0.1;

    let names_tt_2 = names.clone();
    let result = run_optimize_loop(
      &mut graph,
      sparse_partitions,
      dense_partitions,
      max_iter,
      dp,
      damping,
      method,
      false,
      TopologyOps::default(), branch_lengths, &names_tt_2
    )?;
    let sparse_partitions = result.sparse_partitions;
    let dense_partitions = result.dense_partitions;

    assert!(
      result.stopped_at.is_some(),
      "Damped optimization did not stop within {max_iter} iterations"
    );

    let (dense_partitions, sparse_partitions, final_lh) = compute_total_lh(&graph, dense_partitions, sparse_partitions, &result.branch_lengths)?;
    assert!(
      final_lh > -73.0 && final_lh < -72.0,
      "Final log-lh {final_lh:.6} outside expected range (-73.0, -72.0)"
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
  fn test_damped_optimization_does_not_regress(#[case] method: BranchOptMethod) -> Result<(), Report> {
    let aln = simple_alignment()?;
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let (dense_partitions, sparse_partitions) = setup_partitions(&graph, &names, &aln, &mut branch_lengths)?;

    let (dense_partitions, sparse_partitions, initial_lh) = compute_total_lh(&graph, dense_partitions, sparse_partitions, &branch_lengths)?;

    let dp = 0.0;
    let names_tt_1 = names.clone();
    let result = run_optimize_loop(
      &mut graph,
      sparse_partitions,
      dense_partitions,
      10,
      dp,
      0.75,
      method,
      false,
      TopologyOps::default(), branch_lengths, &names_tt_1
    )?;
    let sparse_partitions = result.sparse_partitions;
    let dense_partitions = result.dense_partitions;

    let (dense_partitions, sparse_partitions, final_lh) = compute_total_lh(&graph, dense_partitions, sparse_partitions, &result.branch_lengths)?;
    assert!(
      final_lh >= initial_lh,
      "Damped optimization regressed: {initial_lh:.6} -> {final_lh:.6}"
    );
    Ok(())
  }
}
