#[cfg(test)]
mod tests {
  use super::super::test_convergence_support::tests::{
    TREE_NEWICK, compute_total_lh, setup_reconstruction, simple_alignment,
  };
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::optimize::dispatch::run_optimize_mixed;
  use crate::optimize::gather::{gather_edge_contributions, gather_edge_indel_counts};
  use crate::optimize::params::BranchOptMethod;
  use crate::pretty_assert_ulps_eq;
  use eyre::Report;
  use rstest::rstest;
  use treetime_io::nwk::nwk_read_str;

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
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;

    let mut reconstruction = setup_reconstruction(&graph, &names, &aln, &mut branch_lengths)?;

    let mut lh_history = Vec::with_capacity(20);

    for i in 0..20 {
      let total_length = reconstruction.sequence_length();
      let contributions = gather_edge_contributions(&graph, &reconstruction)?;
      let indel_counts = gather_edge_indel_counts(&graph, &reconstruction);
      run_optimize_mixed(&graph, total_length, &contributions, &indel_counts, method, &mut branch_lengths)?;
      let (reconstruction_updated, log_lh) =
        reconstruction.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
      reconstruction = reconstruction_updated;
      let lh = log_lh.value();

      lh_history.push(lh);

      for edge in graph.get_edges() {
        let branch_length = branch_lengths[&edge.key()];
        if let Some(bl) = branch_length {
          assert!(bl >= 0.0, "Branch length should be non-negative at iter {i}: {bl}");
          assert!(bl < 10.0, "Branch length too large at iter {i}: {bl}");
        }
      }
    }

    let last_5: Vec<f64> = lh_history.iter().rev().take(5).copied().collect();
    let mean = last_5.iter().sum::<f64>() / 5.0;
    let variance = last_5.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / 5.0;
    assert!(
      variance < 1.0,
      "Optimization should stabilize: variance of last 5 iterations = {variance}"
    );

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

    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let graph1_names = nwk_parsed.names();
    let graph1 = nwk_parsed.graph;
    let mut branch_lengths1 = nwk_parsed.branch_lengths;
    let mut reconstruction1 = setup_reconstruction(&graph1, &graph1_names, &aln, &mut branch_lengths1)?;

    for _ in 0..10 {
      let total_length = reconstruction1.sequence_length();
      let contributions = gather_edge_contributions(&graph1, &reconstruction1)?;
      let indel_counts = gather_edge_indel_counts(&graph1, &reconstruction1);
      run_optimize_mixed(&graph1, total_length, &contributions, &indel_counts, method, &mut branch_lengths1)?;
      (reconstruction1, _) = reconstruction1.marginal_update(&graph1, &branch_lengths_or_zero(&branch_lengths1))?;
    }

    let (_, lh1) = compute_total_lh(&graph1, reconstruction1, &branch_lengths1)?;

    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let graph2_names = nwk_parsed.names();
    let graph2 = nwk_parsed.graph;
    let mut branch_lengths2 = nwk_parsed.branch_lengths;
    let mut reconstruction2 = setup_reconstruction(&graph2, &graph2_names, &aln, &mut branch_lengths2)?;

    for _ in 0..10 {
      let total_length = reconstruction2.sequence_length();
      let contributions = gather_edge_contributions(&graph2, &reconstruction2)?;
      let indel_counts = gather_edge_indel_counts(&graph2, &reconstruction2);
      run_optimize_mixed(&graph2, total_length, &contributions, &indel_counts, method, &mut branch_lengths2)?;
      (reconstruction2, _) = reconstruction2.marginal_update(&graph2, &branch_lengths_or_zero(&branch_lengths2))?;
    }

    let (_, lh2) = compute_total_lh(&graph2, reconstruction2, &branch_lengths2)?;

    pretty_assert_ulps_eq!(lh1, lh2, max_ulps = 100);

    Ok(())
  }
}
