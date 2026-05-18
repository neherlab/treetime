#[cfg(test)]
mod tests {
  use crate::ancestral::marginal::update_marginal;
  use crate::optimize::params::BranchOptMethod;
  use crate::optimize::dispatch::run_optimize_mixed;
  use crate::pretty_assert_ulps_eq;
  use crate::partition::payload::ancestral::GraphAncestral;
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
  fn test_optimization_converges_with_valid_branch_lengths(#[case] method: BranchOptMethod) -> Result<(), Report> {
    let aln = simple_alignment()?;
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;

    let (dense_partitions, sparse_partitions, mixed_partitions) = setup_partitions(&graph, &aln)?;

    let mut lh_history = Vec::with_capacity(20);

    // Run optimization iterations
    for i in 0..20 {
      run_optimize_mixed(&graph, &mixed_partitions, method)?;
      let lh = update_marginal(&graph, &dense_partitions)? + update_marginal(&graph, &sparse_partitions)?;

      lh_history.push(lh);

      // After each iteration, all branch lengths should be non-negative and bounded
      for edge in graph.get_edges() {
        let branch_length = edge.read_arc().payload().read_arc().branch_length();
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
    let graph1: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (dense_partitions1, sparse_partitions1, mixed_partitions1) = setup_partitions(&graph1, &aln)?;

    for _ in 0..10 {
      run_optimize_mixed(&graph1, &mixed_partitions1, method)?;
      update_marginal(&graph1, &dense_partitions1)?;
      update_marginal(&graph1, &sparse_partitions1)?;
    }

    let lh1 = compute_total_lh(&graph1, &dense_partitions1, &sparse_partitions1)?;

    // Run optimization on second independent graph
    let graph2: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (dense_partitions2, sparse_partitions2, mixed_partitions2) = setup_partitions(&graph2, &aln)?;

    for _ in 0..10 {
      run_optimize_mixed(&graph2, &mixed_partitions2, method)?;
      update_marginal(&graph2, &dense_partitions2)?;
      update_marginal(&graph2, &sparse_partitions2)?;
    }

    let lh2 = compute_total_lh(&graph2, &dense_partitions2, &sparse_partitions2)?;

    // Both runs should converge to same likelihood
    pretty_assert_ulps_eq!(lh1, lh2, max_ulps = 100);

    Ok(())
  }
}
