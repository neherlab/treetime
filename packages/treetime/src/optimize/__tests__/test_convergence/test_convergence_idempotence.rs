#[cfg(test)]
mod tests {
  use crate::ancestral::marginal::profile_branch_lengths;
  use crate::optimize::dispatch::run_optimize_mixed;
  use crate::optimize::params::BranchOptMethod;
  use crate::optimize::run_loop::{OptimizeReadouts, marginal_update_dense, marginal_update_sparse};
  use crate::pretty_assert_ulps_eq;
  use eyre::Report;
  use rstest::rstest;
  use treetime_graph::graph::Graph;
  use treetime_io::nwk::{NwkParse, nwk_read_str};

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
  fn test_optimization_converges_with_valid_branch_lengths(#[case] method: BranchOptMethod) -> Result<(), Report> {
    let aln = simple_alignment()?;
    let NwkParse { graph, names, mut branch_lengths, .. } = nwk_read_str(TREE_NEWICK)?;
    let graph: Graph = graph;

    let (mut dense_partitions, mut sparse_partitions) = setup_partitions(&graph, &names, &aln, &mut branch_lengths)?;

    let mut lh_history = Vec::with_capacity(20);

    // Run optimization iterations
    for i in 0..20 {
      run_optimize_mixed(&graph, &OptimizeReadouts::new(&dense_partitions, &sparse_partitions).view(), method, &mut branch_lengths)?;
      let (dense_partitions_updated, dense_lh) =
        marginal_update_dense(&graph, &profile_branch_lengths(&branch_lengths), dense_partitions)?;
      let (sparse_partitions_updated, sparse_lh) =
        marginal_update_sparse(&graph, &profile_branch_lengths(&branch_lengths), sparse_partitions)?;
      dense_partitions = dense_partitions_updated;
      sparse_partitions = sparse_partitions_updated;
      let lh = dense_lh.value() + sparse_lh.value();

      lh_history.push(lh);

      // After each iteration, all branch lengths should be non-negative and bounded
      for edge in graph.get_edges() {
        let branch_length = branch_lengths[&edge.read_arc().key()];
        if let Some(bl) = branch_length {
          assert!(bl >= 0.0, "Branch length should be non-negative at iter {i}: {bl}");
          assert!(bl < 10.0, "Branch length too large at iter {i}: {bl}");
        }
      }
    }

    // Check convergence: variance over last 5 iterations should be small
    let last_5: Vec<f64> = lh_history.iter().rev().take(5).copied().collect();
    let mean = last_5.iter().sum::<f64>() / 5.0;
    let variance = last_5.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / 5.0;
    assert!(
      variance < 1.0,
      "Optimization should stabilize: variance of last 5 iterations = {variance}"
    );

    // Final likelihood should be in expected range
    let final_lh = lh_history[19];
    assert!(final_lh < 0.0, "Final log-LH should be negative: {final_lh}");
    assert!(final_lh > -100.0, "Final log-LH unreasonably low: {final_lh}");

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
  fn test_second_optimization_produces_same_likelihood(#[case] method: BranchOptMethod) -> Result<(), Report> {
    let aln = simple_alignment()?;

    // Run optimization on first graph
    let NwkParse { graph: graph1, names: graph1_names, branch_lengths: mut branch_lengths1, .. } = nwk_read_str(TREE_NEWICK)?;
    let (mut dense_partitions1, mut sparse_partitions1) = setup_partitions(&graph1, &graph1_names, &aln, &mut branch_lengths1)?;

    for _ in 0..10 {
      run_optimize_mixed(&graph1, &OptimizeReadouts::new(&dense_partitions1, &sparse_partitions1).view(), method, &mut branch_lengths1)?;
      (dense_partitions1, _) = marginal_update_dense(&graph1, &profile_branch_lengths(&branch_lengths1), dense_partitions1)?;
      (sparse_partitions1, _) = marginal_update_sparse(&graph1, &profile_branch_lengths(&branch_lengths1), sparse_partitions1)?;
    }

    let (dense_partitions1, sparse_partitions1, lh1) = compute_total_lh(&graph1, dense_partitions1, sparse_partitions1, &branch_lengths1)?;

    // Run optimization on second independent graph
    let NwkParse { graph: graph2, names: graph2_names, branch_lengths: mut branch_lengths2, .. } = nwk_read_str(TREE_NEWICK)?;
    let (mut dense_partitions2, mut sparse_partitions2) = setup_partitions(&graph2, &graph2_names, &aln, &mut branch_lengths2)?;

    for _ in 0..10 {
      run_optimize_mixed(&graph2, &OptimizeReadouts::new(&dense_partitions2, &sparse_partitions2).view(), method, &mut branch_lengths2)?;
      (dense_partitions2, _) = marginal_update_dense(&graph2, &profile_branch_lengths(&branch_lengths2), dense_partitions2)?;
      (sparse_partitions2, _) = marginal_update_sparse(&graph2, &profile_branch_lengths(&branch_lengths2), sparse_partitions2)?;
    }

    let (dense_partitions2, sparse_partitions2, lh2) = compute_total_lh(&graph2, dense_partitions2, sparse_partitions2, &branch_lengths2)?;

    // Both runs should converge to same likelihood
    pretty_assert_ulps_eq!(lh1, lh2, max_ulps = 100);

    Ok(())
  }
}
