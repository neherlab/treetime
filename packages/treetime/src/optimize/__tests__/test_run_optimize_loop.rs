#![allow(
  clippy::as_conversions,
  reason = "test and benchmark code: index and expected-value casts, property-style tests over thread_rng inputs (seeding is a separate test-quality follow-up), and scratch collections"
)]

#[cfg(test)]
mod tests {
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::optimize::__tests__::test_convergence::test_convergence_support::tests::{
    TREE_NEWICK, setup_reconstruction, simple_alignment,
  };
  use crate::optimize::params::{BranchOptMethod, TopologyOps};
  use crate::optimize::run_loop::{ConvergenceReason, run_optimize_loop};
  use crate::test_utils::{deletion, sparse_reconstruction_mut};
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use helpers::manual_total_indel_log_lh;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::Seq;

  #[test]
  fn test_run_optimize_loop_records_lh_history() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let mut graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let reconstruction = setup_reconstruction(&graph, &names, &aln, &mut branch_lengths)?;

    let max_iter = 5;
    let names_tt_6 = names.clone();
    let result = run_optimize_loop(
      &mut graph,
      reconstruction,
      max_iter,
      0.0,
      0.75,
      BranchOptMethod::BrentSqrt,
      false,
      TopologyOps::default(),
      branch_lengths,
      &names_tt_6,
    )?;

    assert!(result.lh_history.len() >= 1);
    assert!(result.lh_history.len() <= max_iter);
    Ok(())
  }

  #[test]
  fn test_run_optimize_loop_records_joint_likelihood_with_sparse_indels() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let mut graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut reconstruction = setup_reconstruction(&graph, &names, &aln, &mut branch_lengths)?;

    let first_edge_key = graph.get_edges().collect::<Vec<_>>()[0].key();
    branch_lengths.insert(first_edge_key, Some(0.1));
    sparse_reconstruction_mut(&mut reconstruction)
      .partition
      .obs_edges
      .get_mut(&first_edge_key)
      .unwrap()
      .indels = vec![deletion((0, 3), Seq::try_from_str("ACG")?)];

    let (reconstruction, sparse_lh) =
      reconstruction.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
    let sparse_lh = sparse_lh.value();
    let indel_lh = manual_total_indel_log_lh(&graph, &reconstruction, &branch_lengths);
    let expected_total_lh = sparse_lh + indel_lh;

    let names_tt_5 = names.clone();
    let result = run_optimize_loop(
      &mut graph,
      reconstruction,
      1,
      0.0,
      0.75,
      BranchOptMethod::BrentSqrt,
      false,
      TopologyOps::default(),
      branch_lengths,
      &names_tt_5,
    )?;

    assert_eq!(1, result.lh_history.len());
    assert_abs_diff_eq!(result.lh_history[0].value(), expected_total_lh, epsilon = 1e-10);
    Ok(())
  }

  #[test]
  fn test_run_optimize_loop_breaks_on_convergence() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let mut graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let reconstruction = setup_reconstruction(&graph, &names, &aln, &mut branch_lengths)?;

    let max_iter = 50;
    let dp = f64::INFINITY;

    let names_tt_4 = names.clone();
    let result = run_optimize_loop(
      &mut graph,
      reconstruction,
      max_iter,
      dp,
      0.0,
      BranchOptMethod::BrentSqrt,
      false,
      TopologyOps::default(),
      branch_lengths,
      &names_tt_4,
    )?;

    assert_eq!(Some((1, ConvergenceReason::Converged)), result.stopped_at);
    assert_eq!(2, result.lh_history.len());
    Ok(())
  }

  #[test]
  fn test_run_optimize_loop_zero_max_iter_is_noop() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let mut graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let reconstruction = setup_reconstruction(&graph, &names, &aln, &mut branch_lengths)?;

    let names_tt_3 = names.clone();
    let result = run_optimize_loop(
      &mut graph,
      reconstruction,
      0,
      1e-2,
      0.75,
      BranchOptMethod::BrentSqrt,
      false,
      TopologyOps::default(),
      branch_lengths,
      &names_tt_3,
    )?;

    assert!(result.lh_history.is_empty());
    assert!(result.stopped_at.is_none());
    Ok(())
  }

  #[test]
  fn test_run_optimize_loop_all_likelihoods_finite() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let mut graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let reconstruction = setup_reconstruction(&graph, &names, &aln, &mut branch_lengths)?;

    let names_tt_2 = names.clone();
    let result = run_optimize_loop(
      &mut graph,
      reconstruction,
      10,
      0.0,
      0.75,
      BranchOptMethod::BrentSqrt,
      false,
      TopologyOps::default(),
      branch_lengths,
      &names_tt_2,
    )?;

    for (i, lh) in result.lh_history.iter().map(|log_lh| log_lh.value()).enumerate() {
      assert!(lh.is_finite(), "Iteration {i}: log-likelihood must be finite, got {lh}");
    }
    Ok(())
  }

  #[test]
  fn test_run_optimize_loop_improves_likelihood() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let mut graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut reconstruction = setup_reconstruction(&graph, &names, &aln, &mut branch_lengths)?;
    let first_edge_key = graph.get_edges().collect::<Vec<_>>()[0].key();
    branch_lengths.insert(first_edge_key, Some(0.1));
    sparse_reconstruction_mut(&mut reconstruction)
      .partition
      .obs_edges
      .get_mut(&first_edge_key)
      .unwrap()
      .indels = vec![deletion((0, 3), Seq::try_from_str("ACG")?)];

    let (reconstruction, initial_sparse_lh) =
      reconstruction.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
    let initial_lh = initial_sparse_lh.value() + manual_total_indel_log_lh(&graph, &reconstruction, &branch_lengths);

    let names_tt_1 = names.clone();
    let result = run_optimize_loop(
      &mut graph,
      reconstruction,
      10,
      0.0,
      0.75,
      BranchOptMethod::BrentSqrt,
      false,
      TopologyOps::default(),
      branch_lengths,
      &names_tt_1,
    )?;
    let branch_lengths = result.branch_lengths;

    let (reconstruction, final_sparse_lh) = result
      .reconstruction
      .marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
    let final_lh = final_sparse_lh.value() + manual_total_indel_log_lh(&graph, &reconstruction, &branch_lengths);

    assert!(
      final_lh >= initial_lh,
      "Loop regressed likelihood: {initial_lh:.6} -> {final_lh:.6}"
    );
    Ok(())
  }

  mod helpers {
    use crate::partition::marginal::reconstruction::MarginalReconstruction;
    use crate::test_utils::sparse_reconstruction;
    use statrs::function::factorial::ln_factorial;
    use std::collections::BTreeMap;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::graph::Graph;

    fn manual_indel_count_on_edge(reconstruction: &MarginalReconstruction, edge_key: GraphEdgeKey) -> usize {
      sparse_reconstruction(reconstruction)
        .partition
        .obs_edges
        .get(&edge_key)
        .map_or(0, |edge| edge.indels.len())
    }

    fn manual_poisson_indel_log_lh(k: usize, mu: f64, t: f64) -> f64 {
      if k > 0 && t <= 0.0 {
        return f64::NEG_INFINITY;
      }
      if mu == 0.0 {
        return if k == 0 { 0.0 } else { f64::NEG_INFINITY };
      }
      if k == 0 {
        return -mu * t;
      }

      let lambda = mu * t;
      (k as f64) * lambda.ln() - lambda - ln_factorial(k as u64)
    }

    pub(super) fn manual_total_indel_log_lh(
      graph: &Graph,
      reconstruction: &MarginalReconstruction,
      branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
    ) -> f64 {
      let total_indels: usize = graph
        .get_edges()
        .map(|edge_ref| manual_indel_count_on_edge(reconstruction, edge_ref.key()))
        .sum();
      let total_branch_length: f64 = graph
        .get_edges()
        .map(|edge_ref| branch_lengths.get(&edge_ref.key()).copied().flatten().unwrap_or(0.0))
        .sum();
      let indel_rate = if total_indels > 0 && total_branch_length > 0.0 {
        total_indels as f64 / total_branch_length
      } else {
        0.0
      };

      graph
        .get_edges()
        .map(|edge_ref| {
          let edge_key = edge_ref.key();
          let branch_length = branch_lengths.get(&edge_key).copied().flatten().unwrap_or(0.0);
          let indel_count = manual_indel_count_on_edge(reconstruction, edge_key);
          manual_poisson_indel_log_lh(indel_count, indel_rate, branch_length)
        })
        .sum()
    }
  }
}
