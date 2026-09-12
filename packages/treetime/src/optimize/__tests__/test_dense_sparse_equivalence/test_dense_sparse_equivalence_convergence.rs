#[cfg(test)]
mod tests {
  use crate::ancestral::marginal::{marginal_update, profile_branch_lengths};
  use crate::optimize::dispatch::run_optimize_mixed;
  use crate::optimize::params::BranchOptMethod;
  use crate::optimize::run_loop::optimize_partition_view;
  use eyre::Report;
  use rstest::rstest;
  use treetime_graph::graph::Graph;
  use treetime_io::nwk::{NwkParse, nwk_read_str};

  use super::super::test_dense_sparse_equivalence_support::tests::{
    TREE_NEWICK, gap_free_alignment, setup_dense_only, setup_sparse_only,
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
  fn test_dense_optimization_converges(#[case] method: BranchOptMethod) -> Result<(), Report> {
    let aln = gap_free_alignment()?;
    let NwkParse { graph, names, mut branch_lengths, .. } = nwk_read_str(TREE_NEWICK)?;
    let graph: Graph = graph;
    let mut partitions = setup_dense_only(&graph, &names, &aln, &branch_lengths)?;

    let initial_lh = marginal_update(&graph, &profile_branch_lengths(&branch_lengths), &mut partitions)?.value();
    let mut lh_history = vec![initial_lh];

    for _ in 0..50 {
      run_optimize_mixed(&graph, &optimize_partition_view(&partitions, &[]), method, &mut branch_lengths)?;
      let lh = marginal_update(&graph, &profile_branch_lengths(&branch_lengths), &mut partitions)?.value();
      lh_history.push(lh);
    }

    let final_lh = match lh_history.last() {
      Some(final_lh) => *final_lh,
      None => unreachable!("likelihood history always contains the initial value"),
    };

    // Final log-LH should be in expected range
    assert!(
      final_lh > -100.0 && final_lh < -10.0,
      "Final log-LH {final_lh} should be in range [-100, -10]"
    );

    // Optimization should improve or maintain likelihood overall
    assert!(
      final_lh >= initial_lh - 1.0,
      "Optimization should not significantly decrease likelihood: initial={initial_lh}, final={final_lh}"
    );

    // Check convergence: variance of last 5 iterations should be small
    let last_5: Vec<f64> = lh_history.iter().rev().take(5).copied().collect();
    let mean: f64 = last_5.iter().sum::<f64>() / last_5.len() as f64;
    let variance: f64 = last_5.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / last_5.len() as f64;
    assert!(
      variance < 1.0,
      "Optimization should stabilize: variance of last 5 iterations = {variance}"
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
  fn test_sparse_optimization_converges(#[case] method: BranchOptMethod) -> Result<(), Report> {
    let aln = gap_free_alignment()?;
    let NwkParse { graph, names, mut branch_lengths, .. } = nwk_read_str(TREE_NEWICK)?;
    let graph: Graph = graph;
    let mut partitions = setup_sparse_only(&graph, &names, &aln, &branch_lengths)?;

    let initial_lh = marginal_update(&graph, &profile_branch_lengths(&branch_lengths), &mut partitions)?.value();
    let mut lh_history = vec![initial_lh];

    for _ in 0..50 {
      run_optimize_mixed(&graph, &optimize_partition_view(&[], &partitions), method, &mut branch_lengths)?;
      let lh = marginal_update(&graph, &profile_branch_lengths(&branch_lengths), &mut partitions)?.value();
      lh_history.push(lh);
    }

    let final_lh = match lh_history.last() {
      Some(final_lh) => *final_lh,
      None => unreachable!("likelihood history always contains the initial value"),
    };

    // Final log-LH should be in expected range
    assert!(
      final_lh > -100.0 && final_lh < -10.0,
      "Final log-LH {final_lh} should be in range [-100, -10]"
    );

    // All iterations should stay in valid range
    assert!(
      lh_history.iter().all(|lh| *lh > -200.0 && *lh < 0.0),
      "All log-LH values should be in valid range [-200, 0]: {lh_history:?}"
    );

    // Best likelihood achieved should be reasonable
    let best_lh = lh_history.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    assert!(best_lh > -50.0, "Best log-LH {best_lh} should be better than -50");

    Ok(())
  }
}
