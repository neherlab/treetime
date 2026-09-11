#[cfg(test)]
mod tests {
  use crate::ancestral::marginal::{marginal_update, profile_branch_lengths};
  use crate::optimize::dispatch::run_optimize_mixed;
  use crate::optimize::params::BranchOptMethod;
  use crate::payload::ancestral::GraphAncestral;
  use crate::pretty_assert_ulps_eq;
  use eyre::Report;
  use rstest::rstest;
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
    let graph: GraphAncestral = graph;

    let (dense_partitions, sparse_partitions, mixed_partitions) = setup_partitions(&graph, &names, &aln, &mut branch_lengths)?;

    let mut lh_history = Vec::with_capacity(20);

    // Run optimization iterations
    for i in 0..20 {
      run_optimize_mixed(&graph, &mixed_partitions, method, &mut branch_lengths)?;
      let lh = marginal_update(&graph, &profile_branch_lengths(&branch_lengths), &dense_partitions)?.value() + marginal_update(&graph, &profile_branch_lengths(&branch_lengths), &sparse_partitions)?.value();

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
    let (dense_partitions1, sparse_partitions1, mixed_partitions1) = setup_partitions(&graph1, &graph1_names, &aln, &mut branch_lengths1)?;

    for _ in 0..10 {
      run_optimize_mixed(&graph1, &mixed_partitions1, method, &mut branch_lengths1)?;
      marginal_update(&graph1, &profile_branch_lengths(&branch_lengths1), &dense_partitions1)?.value();
      marginal_update(&graph1, &profile_branch_lengths(&branch_lengths1), &sparse_partitions1)?.value();
    }

    let lh1 = compute_total_lh(&graph1, &dense_partitions1, &sparse_partitions1, &branch_lengths1)?;

    // Run optimization on second independent graph
    let NwkParse { graph: graph2, names: graph2_names, branch_lengths: mut branch_lengths2, .. } = nwk_read_str(TREE_NEWICK)?;
    let (dense_partitions2, sparse_partitions2, mixed_partitions2) = setup_partitions(&graph2, &graph2_names, &aln, &mut branch_lengths2)?;

    for _ in 0..10 {
      run_optimize_mixed(&graph2, &mixed_partitions2, method, &mut branch_lengths2)?;
      marginal_update(&graph2, &profile_branch_lengths(&branch_lengths2), &dense_partitions2)?.value();
      marginal_update(&graph2, &profile_branch_lengths(&branch_lengths2), &sparse_partitions2)?.value();
    }

    let lh2 = compute_total_lh(&graph2, &dense_partitions2, &sparse_partitions2, &branch_lengths2)?;

    // Both runs should converge to same likelihood
    pretty_assert_ulps_eq!(lh1, lh2, max_ulps = 100);

    Ok(())
  }
}
