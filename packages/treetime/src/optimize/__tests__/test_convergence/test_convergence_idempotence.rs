#[cfg(test)]
mod tests {
  use super::super::test_convergence_support::tests::{
    TREE_NEWICK, compute_total_lh, setup_reconstruction, simple_alignment,
  };
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::optimize::dispatch::optimize_edge_lengths;
  use crate::optimize::gather::gather_edge_contributions;
  use crate::optimize::indel::IndelPolicy;
  use crate::optimize::params::BranchOptMethod;
  use crate::pretty_assert_ulps_eq;
  use eyre::Report;
  use rstest::rstest;
  use treetime_io::nwk::nwk_read;

  #[rustfmt::skip]
  #[rstest]
  #[case::newton(BranchOptMethod::Newton)]
  #[case::brent( BranchOptMethod::Brent)]
  #[trace]
  fn test_optimization_converges_with_valid_branch_lengths(#[case] method: BranchOptMethod) -> Result<(), Report> {
    let aln = simple_alignment()?;
    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let names = nwk_parsed.names.clone();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;

    let mut reconstruction = setup_reconstruction(&graph, &names, &aln, &mut branch_lengths)?;

    let mut lh_history = Vec::with_capacity(20);

    for i in 0..20 {
      let total_length = reconstruction.total_length();
      let contributions = gather_edge_contributions(&graph, &reconstruction)?;
      let indels = IndelPolicy::Include.edge_indels(&graph, &reconstruction, &branch_lengths)?;
      optimize_edge_lengths(&graph, total_length, &contributions, &indels, method, &mut branch_lengths)?;
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
  #[case::newton(BranchOptMethod::Newton)]
  #[case::brent( BranchOptMethod::Brent)]
  #[trace]
  fn test_second_optimization_produces_same_likelihood(#[case] method: BranchOptMethod) -> Result<(), Report> {
    let aln = simple_alignment()?;

    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let graph1_names = nwk_parsed.names.clone();
    let graph1 = nwk_parsed.graph;
    let mut branch_lengths1 = nwk_parsed.branch_lengths;
    let mut reconstruction1 = setup_reconstruction(&graph1, &graph1_names, &aln, &mut branch_lengths1)?;

    for _ in 0..10 {
      let total_length = reconstruction1.total_length();
      let contributions = gather_edge_contributions(&graph1, &reconstruction1)?;
      let indels = IndelPolicy::Include.edge_indels(&graph1, &reconstruction1, &branch_lengths1)?;
      optimize_edge_lengths(&graph1, total_length, &contributions, &indels, method, &mut branch_lengths1)?;
      (reconstruction1, _) = reconstruction1.marginal_update(&graph1, &branch_lengths_or_zero(&branch_lengths1))?;
    }

    let (_, lh1) = compute_total_lh(&graph1, reconstruction1, &branch_lengths1)?;

    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let graph2_names = nwk_parsed.names.clone();
    let graph2 = nwk_parsed.graph;
    let mut branch_lengths2 = nwk_parsed.branch_lengths;
    let mut reconstruction2 = setup_reconstruction(&graph2, &graph2_names, &aln, &mut branch_lengths2)?;

    for _ in 0..10 {
      let total_length = reconstruction2.total_length();
      let contributions = gather_edge_contributions(&graph2, &reconstruction2)?;
      let indels = IndelPolicy::Include.edge_indels(&graph2, &reconstruction2, &branch_lengths2)?;
      optimize_edge_lengths(&graph2, total_length, &contributions, &indels, method, &mut branch_lengths2)?;
      (reconstruction2, _) = reconstruction2.marginal_update(&graph2, &branch_lengths_or_zero(&branch_lengths2))?;
    }

    let (_, lh2) = compute_total_lh(&graph2, reconstruction2, &branch_lengths2)?;

    pretty_assert_ulps_eq!(lh1, lh2, max_ulps = 100);

    Ok(())
  }
}
