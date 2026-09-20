#[cfg(test)]
mod tests {
  use crate::ancestral::marginal::branch_lengths_or_zero;
  use crate::optimize::__tests__::test_convergence::test_convergence_support::tests::{
    TREE_NEWICK, setup_partitions, simple_alignment,
  };
  use crate::optimize::iteration::{DAMPING_FLOOR, apply_damping};
  use crate::optimize::params::{BranchOptMethod, TopologyOps};
  use crate::optimize::run_loop::{
    ConvergenceReason, marginal_update_dense, marginal_update_sparse, run_optimize_loop,
  };
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use num_traits::pow::pow;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_io::nwk::nwk_read_str;

  #[rustfmt::skip]
  #[rstest]
  #[case::iter_100(  100, DAMPING_FLOOR)]
  #[case::iter_500(  500, DAMPING_FLOOR)]
  #[case::iter_1000(1000, DAMPING_FLOOR)]
  #[trace]
  fn test_convergence_conditions_damping_floor_at_high_iteration(#[case] iteration: usize,
    #[case] expected_old_weight: f64,
  ) -> Result<(), Report> {
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
  fn test_convergence_conditions_damping_uses_exponential_below_crossover() -> Result<(), Report> {
    let damping = 0.75;
    let iteration = 5;
    let expected_old_weight = pow(damping, iteration + 1);
    assert!(expected_old_weight > DAMPING_FLOOR);

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
  fn test_convergence_conditions_converged_reason() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let (dense_partitions, sparse_partitions) = setup_partitions(&graph, &names, &aln, &mut branch_lengths)?;

    let names_tt_6 = names.clone();
    let result = run_optimize_loop(
      &mut graph,
      sparse_partitions,
      dense_partitions,
      20,
      0.1,
      0.75,
      BranchOptMethod::BrentSqrt,
      false,
      TopologyOps::default(),
      branch_lengths,
      &names_tt_6,
    )?;
    let sparse_partitions = result.sparse_partitions;
    let dense_partitions = result.dense_partitions;

    let (iter, reason) = result.stopped_at.expect("loop should have stopped");
    assert!(
      reason == ConvergenceReason::Converged || reason == ConvergenceReason::Oscillating,
      "Expected Converged or Oscillating on toy tree, got {reason:?} at iteration {iter}"
    );
    Ok(())
  }

  #[test]
  fn test_convergence_conditions_worsened_reverts_to_best() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let (dense_partitions, sparse_partitions) = setup_partitions(&graph, &names, &aln, &mut branch_lengths)?;

    let names_tt_5 = names.clone();
    let result = run_optimize_loop(
      &mut graph,
      sparse_partitions,
      dense_partitions,
      50,
      0.0,
      0.0,
      BranchOptMethod::BrentSqrt,
      false,
      TopologyOps::default(),
      branch_lengths,
      &names_tt_5,
    )?;
    let sparse_partitions = result.sparse_partitions;
    let dense_partitions = result.dense_partitions;

    match result.stopped_at {
      Some((iter, ConvergenceReason::Worsened)) => {
        assert!(iter >= 2, "Worsened should not fire before iteration 2, got {iter}");
        let best_lh = result
          .lh_history
          .iter()
          .map(|log_lh| log_lh.value())
          .fold(f64::NEG_INFINITY, f64::max);
        let trigger_lh = result.lh_history[iter].value();
        assert!(
          trigger_lh < best_lh,
          "Trigger LH ({trigger_lh:.6}) should be less than best LH ({best_lh:.6})"
        );
      },
      other => panic!("Expected Worsened on undamped toy tree, got {other:?}"),
    }
    Ok(())
  }

  #[test]
  fn test_convergence_conditions_worsened_rollback_reproduces_best_lh() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let (dense_partitions, sparse_partitions) = setup_partitions(&graph, &names, &aln, &mut branch_lengths)?;

    let names_tt_4 = names.clone();
    let result = run_optimize_loop(
      &mut graph,
      sparse_partitions,
      dense_partitions,
      50,
      0.0,
      0.0,
      BranchOptMethod::BrentSqrt,
      true,
      TopologyOps::default(),
      branch_lengths,
      &names_tt_4,
    )?;
    let sparse_partitions = result.sparse_partitions;
    let dense_partitions = result.dense_partitions;

    let (_iter, reason) = result.stopped_at.expect("loop should stop");
    assert_eq!(ConvergenceReason::Worsened, reason);

    let best_lh = result
      .lh_history
      .iter()
      .map(|log_lh| log_lh.value())
      .fold(f64::NEG_INFINITY, f64::max);

    let marginal_bl = branch_lengths_or_zero(&result.branch_lengths);
    let (sparse_partitions, sparse_lh) = marginal_update_sparse(&graph, &marginal_bl, sparse_partitions)?;
    let sparse_lh = sparse_lh.value();
    let (dense_partitions, dense_lh) = marginal_update_dense(&graph, &marginal_bl, dense_partitions)?;
    let dense_lh = dense_lh.value();
    assert_abs_diff_eq!(sparse_lh + dense_lh, best_lh, epsilon = 1e-9);
    Ok(())
  }

  #[test]
  fn test_convergence_conditions_oscillation_detection() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let (dense_partitions, sparse_partitions) = setup_partitions(&graph, &names, &aln, &mut branch_lengths)?;

    let names_tt_3 = names.clone();
    let result = run_optimize_loop(
      &mut graph,
      sparse_partitions,
      dense_partitions,
      50,
      1.0,
      0.75,
      BranchOptMethod::BrentSqrt,
      false,
      TopologyOps::default(),
      branch_lengths,
      &names_tt_3,
    )?;
    let sparse_partitions = result.sparse_partitions;
    let dense_partitions = result.dense_partitions;

    let (iter, reason) = result.stopped_at.expect("loop should have stopped");
    assert!(
      reason == ConvergenceReason::Converged || reason == ConvergenceReason::Oscillating,
      "Expected early stop with dp=1.0, got {reason:?} at iteration {iter}"
    );
    assert!(iter < 10, "Expected early stop, got iteration {iter}");
    Ok(())
  }

  #[test]
  fn test_convergence_conditions_exhausts_max_iter() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let (dense_partitions, sparse_partitions) = setup_partitions(&graph, &names, &aln, &mut branch_lengths)?;

    let names_tt_2 = names.clone();
    let result = run_optimize_loop(
      &mut graph,
      sparse_partitions,
      dense_partitions,
      2,
      0.0,
      0.75,
      BranchOptMethod::BrentSqrt,
      false,
      TopologyOps::default(),
      branch_lengths,
      &names_tt_2,
    )?;
    let sparse_partitions = result.sparse_partitions;
    let dense_partitions = result.dense_partitions;

    assert_eq!(2, result.lh_history.len());
    assert!(
      result.stopped_at.is_none(),
      "Should exhaust max_iter=2 without stopping"
    );
    Ok(())
  }

  #[test]
  fn test_convergence_conditions_dense_only_converges() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;

    let (dense_partitions, _sparse_partitions) = setup_partitions(&graph, &names, &aln, &mut branch_lengths)?;

    let empty_sparse = vec![];

    let names_tt_1 = names.clone();
    let result = run_optimize_loop(
      &mut graph,
      empty_sparse,
      dense_partitions,
      10,
      0.1,
      0.75,
      BranchOptMethod::BrentSqrt,
      false,
      TopologyOps::default(),
      branch_lengths,
      &names_tt_1,
    )?;
    let sparse_partitions = result.sparse_partitions;
    let dense_partitions = result.dense_partitions;

    assert!(
      result.stopped_at.is_some(),
      "Dense-only optimization should converge within 10 iterations"
    );
    Ok(())
  }
}
