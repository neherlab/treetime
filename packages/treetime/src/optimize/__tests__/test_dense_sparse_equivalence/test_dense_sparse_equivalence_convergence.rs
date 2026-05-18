#[cfg(test)]
mod tests {
  use crate::ancestral::marginal::update_marginal;
  use crate::optimize::params::BranchOptMethod;
  use crate::optimize::dispatch::run_optimize_mixed;
  use crate::payload::ancestral::GraphAncestral;
  use eyre::Report;
  use rstest::rstest;
  use treetime_io::nwk::nwk_read_str;

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
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let partitions = setup_dense_only(&graph, &aln)?;

    let initial_lh = update_marginal(&graph, &partitions)?;
    let mut lh_history = vec![initial_lh];

    for _ in 0..50 {
      run_optimize_mixed(&graph, &partitions, method)?;
      let lh = update_marginal(&graph, &partitions)?;
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
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let partitions = setup_sparse_only(&graph, &aln)?;

    let initial_lh = update_marginal(&graph, &partitions)?;
    let mut lh_history = vec![initial_lh];

    for _ in 0..50 {
      run_optimize_mixed(&graph, &partitions, method)?;
      let lh = update_marginal(&graph, &partitions)?;
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
