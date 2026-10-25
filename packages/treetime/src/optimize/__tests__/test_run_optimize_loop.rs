#[cfg(test)]
mod tests {
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::optimize::__tests__::test_convergence::test_convergence_support::tests::{
    TREE_NEWICK, setup_dense_reconstruction, setup_reconstruction, simple_alignment,
  };
  use crate::optimize::indel::IndelPolicy;
  use crate::optimize::params::OptimizeLoopSettings;
  use crate::optimize::params::{BranchOptMethod, TopologyOps};
  use crate::optimize::partition_rates::PartitionRates;
  use crate::optimize::run_loop::{ConvergenceReason, run_optimize_loop};
  use crate::progress::{LogLevel, NoopProgress};
  use crate::test_utils::{RecordingLog, SingleReconstruction};
  use crate::test_utils::{deletion, sparse_reconstruction_mut};
  use eyre::Report;
  use helpers::{chain_alignment, chain_alignment_jc69_distance, manual_total_indel_log_lh};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_graph::key_map::EdgeMap;

  use treetime_io::nwk::nwk_read;
  use treetime_primitives::Seq;
  use treetime_utils::pretty_assert_abs_diff_eq;

  #[test]
  fn test_run_optimize_loop_records_lh_history() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let names = nwk_parsed.names.clone();
    let mut graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let reconstruction = setup_reconstruction(&graph, &names, &aln, &mut branch_lengths)?;

    let max_iter = 5;
    let names_tt_6 = names.clone();
    let result = run_optimize_loop(
      &mut graph,
      reconstruction,
      &OptimizeLoopSettings {
        max_iter,
        dp: 0.0,
        damping: 0.75,
        opt_method: BranchOptMethod::Brent,
        indel_policy: IndelPolicy::Include,
        topology_ops: TopologyOps::default(),
        partition_rates: PartitionRates::Estimate,
      },
      branch_lengths,
      names_tt_6,
      &NoopProgress,
    )?;

    assert!(result.lh_history.len() >= 1);
    assert!(result.lh_history.len() <= max_iter + 1);
    Ok(())
  }

  #[test]
  fn test_run_optimize_loop_records_joint_likelihood_with_sparse_indels() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let names = nwk_parsed.names.clone();
    let mut graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut reconstruction = setup_reconstruction(&graph, &names, &aln, &mut branch_lengths)?;

    let first_edge_key = graph.get_edges().collect::<Vec<_>>()[0].key();
    branch_lengths.insert(first_edge_key, Some(0.1));
    sparse_reconstruction_mut(&mut reconstruction)
      .partition
      .obs_edges
      .get_mut(first_edge_key)
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
      &OptimizeLoopSettings {
        max_iter: 1,
        dp: 0.0,
        damping: 0.75,
        opt_method: BranchOptMethod::Brent,
        indel_policy: IndelPolicy::Include,
        topology_ops: TopologyOps::default(),
        partition_rates: PartitionRates::Estimate,
      },
      branch_lengths,
      names_tt_5,
      &NoopProgress,
    )?;

    assert_eq!(2, result.lh_history.len());
    pretty_assert_abs_diff_eq!(result.lh_history[0].value(), expected_total_lh, epsilon = 1e-10);
    Ok(())
  }

  #[test]
  fn test_run_optimize_loop_breaks_on_convergence() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let names = nwk_parsed.names.clone();
    let mut graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let reconstruction = setup_reconstruction(&graph, &names, &aln, &mut branch_lengths)?;

    let max_iter = 50;
    let dp = f64::INFINITY;

    let names_tt_4 = names.clone();
    let result = run_optimize_loop(
      &mut graph,
      reconstruction,
      &OptimizeLoopSettings {
        max_iter,
        dp,
        damping: 0.0,
        opt_method: BranchOptMethod::Brent,
        indel_policy: IndelPolicy::Include,
        topology_ops: TopologyOps::default(),
        partition_rates: PartitionRates::Estimate,
      },
      branch_lengths,
      names_tt_4,
      &NoopProgress,
    )?;

    assert_eq!(Some((1, ConvergenceReason::Converged)), result.stopped_at);
    assert_eq!(2, result.lh_history.len());
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::converged(      (50, f64::INFINITY), "Optimization loop stopped at iteration 2: the log-likelihood changed by less than inf")]
  #[case::iteration_limit((1,  0.0),           "Optimization loop reached the iteration limit of 1")]
  #[case::zero_max_iter(  (0,  1e-2),          "Optimization loop skipped: the iteration limit is 0")]
  #[trace]
  fn test_run_optimize_loop_logs_stop_reason(
    #[case] (max_iter, dp): (usize, f64),
    #[case] expected: &str,
  ) -> Result<(), Report> {
    let aln = simple_alignment()?;
    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let names = nwk_parsed.names.clone();
    let mut graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let reconstruction = setup_reconstruction(&graph, &names, &aln, &mut branch_lengths)?;
    let log = RecordingLog::default();

    run_optimize_loop(
      &mut graph,
      reconstruction,
      &OptimizeLoopSettings {
        max_iter,
        dp,
        damping: 0.0,
        opt_method: BranchOptMethod::Brent,
        indel_policy: IndelPolicy::Include,
        topology_ops: TopologyOps::default(),
        partition_rates: PartitionRates::Estimate,
      },
      branch_lengths,
      names,
      &log,
    )?;

    let stop_messages = log
      .messages(LogLevel::Info)
      .into_iter()
      .filter(|message| message.starts_with("Optimization loop"))
      .collect::<Vec<_>>();
    assert_eq!(vec![expected.to_owned()], stop_messages);
    Ok(())
  }

  #[test]
  fn test_run_optimize_loop_zero_max_iter_is_noop() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let names = nwk_parsed.names.clone();
    let mut graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let reconstruction = setup_reconstruction(&graph, &names, &aln, &mut branch_lengths)?;

    let names_tt_3 = names.clone();
    let result = run_optimize_loop(
      &mut graph,
      reconstruction,
      &OptimizeLoopSettings {
        max_iter: 0,
        dp: 1e-2,
        damping: 0.75,
        opt_method: BranchOptMethod::Brent,
        indel_policy: IndelPolicy::Include,
        topology_ops: TopologyOps::default(),
        partition_rates: PartitionRates::Estimate,
      },
      branch_lengths,
      names_tt_3,
      &NoopProgress,
    )?;

    assert!(result.lh_history.is_empty());
    assert!(result.stopped_at.is_none());
    Ok(())
  }

  #[test]
  fn test_run_optimize_loop_all_likelihoods_finite() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let names = nwk_parsed.names.clone();
    let mut graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let reconstruction = setup_reconstruction(&graph, &names, &aln, &mut branch_lengths)?;

    let names_tt_2 = names.clone();
    let result = run_optimize_loop(
      &mut graph,
      reconstruction,
      &OptimizeLoopSettings {
        max_iter: 10,
        dp: 0.0,
        damping: 0.75,
        opt_method: BranchOptMethod::Brent,
        indel_policy: IndelPolicy::Include,
        topology_ops: TopologyOps::default(),
        partition_rates: PartitionRates::Estimate,
      },
      branch_lengths,
      names_tt_2,
      &NoopProgress,
    )?;

    for (i, lh) in result.lh_history.iter().map(|log_lh| log_lh.value()).enumerate() {
      assert!(lh.is_finite(), "Iteration {i}: log-likelihood must be finite, got {lh}");
    }
    Ok(())
  }

  #[test]
  fn test_run_optimize_loop_improves_likelihood() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let names = nwk_parsed.names.clone();
    let mut graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut reconstruction = setup_reconstruction(&graph, &names, &aln, &mut branch_lengths)?;
    let first_edge_key = graph.get_edges().collect::<Vec<_>>()[0].key();
    branch_lengths.insert(first_edge_key, Some(0.1));
    sparse_reconstruction_mut(&mut reconstruction)
      .partition
      .obs_edges
      .get_mut(first_edge_key)
      .unwrap()
      .indels = vec![deletion((0, 3), Seq::try_from_str("ACG")?)];

    let (reconstruction, initial_sparse_lh) =
      reconstruction.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
    let initial_lh = initial_sparse_lh.value() + manual_total_indel_log_lh(&graph, &reconstruction, &branch_lengths);

    let names_tt_1 = names.clone();
    let result = run_optimize_loop(
      &mut graph,
      reconstruction,
      &OptimizeLoopSettings {
        max_iter: 10,
        dp: 0.0,
        damping: 0.75,
        opt_method: BranchOptMethod::Brent,
        indel_policy: IndelPolicy::Include,
        topology_ops: TopologyOps::default(),
        partition_rates: PartitionRates::Estimate,
      },
      branch_lengths,
      names_tt_1,
      &NoopProgress,
    )?;
    let branch_lengths = result.branch_lengths;

    let (reconstruction, final_sparse_lh) = result
      .partitions
      .marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
    let final_lh = final_sparse_lh.value() + manual_total_indel_log_lh(&graph, &reconstruction, &branch_lengths);

    assert!(
      final_lh >= initial_lh,
      "Loop regressed likelihood: {initial_lh:.6} -> {final_lh:.6}"
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::newton(BranchOptMethod::Newton)]
  #[case::brent( BranchOptMethod::Brent)]
  #[trace]
  fn test_run_optimize_loop_returns_first_iteration_when_likelihood_declines_from_start(
    #[case] method: BranchOptMethod,
  ) -> Result<(), Report> {
    let aln = chain_alignment()?;
    let nwk_parsed = nwk_read(b"((A:0.1,X:0.1)n1:0.1,(B:0.1,Y:0.1)n2:0.1)root;".as_slice())?;
    let names = nwk_parsed.names.clone();
    let mut graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut reconstruction = setup_dense_reconstruction(&graph, &names, &aln, &mut branch_lengths)?;
    reconstruction.gtr_mut().set_mu(1.0)?;
    let distance_above_optimum = 0.003;
    let start_branch_length = (chain_alignment_jc69_distance() + distance_above_optimum) / 4.0;
    let start_branch_lengths: EdgeMap<Option<f64>> = branch_lengths
      .keys()
      .map(|edge_key| (edge_key, Some(start_branch_length)))
      .collect();

    let max_iter = 4;
    let dp = 0.0;
    let damping = 0.0;
    let indel_policy = IndelPolicy::Exclude;

    let result = run_optimize_loop(
      &mut graph,
      reconstruction,
      &OptimizeLoopSettings {
        max_iter,
        dp,
        damping,
        opt_method: method,
        indel_policy,
        topology_ops: TopologyOps::default(),
        partition_rates: PartitionRates::Estimate,
      },
      start_branch_lengths.clone(),
      names.clone(),
      &NoopProgress,
    )?;

    assert!(
      result.lh_history.is_sorted_by(|earlier, later| earlier > later),
      "Likelihood must decline at every iteration for this test to apply: {:?}",
      result.lh_history
    );
    assert_eq!(start_branch_lengths, result.branch_lengths);
    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::Alphabet;
    use crate::partition::set::PartitionSet;
    use crate::test_utils::sparse_reconstruction;
    use eyre::Report;
    use indoc::indoc;
    use statrs::function::factorial::ln_factorial;
    use treetime_graph::key_map::EdgeMap;

    use std::sync::LazyLock;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::graph::Graph;
    use treetime_io::fasta::fasta_read;
    use treetime_primitives::AlignmentRecord;
    use treetime_utils::cast::ToF64Exact;

    static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

    pub(super) fn chain_alignment() -> Result<Vec<AlignmentRecord>, Report> {
      Ok(
        fasta_read(
          indoc! {r#"
            >A
            ACGTACGTACGTACGTACGT
            >B
            ACGTACGTACGTACGTACAA
            >X
            NNNNNNNNNNNNNNNNNNNN
            >Y
            NNNNNNNNNNNNNNNNNNNN
          "#}
          .as_bytes(),
          &*NUC_ALPHABET,
        )?
        .into_iter()
        .map(AlignmentRecord::from)
        .collect(),
      )
    }

    pub(super) fn chain_alignment_jc69_distance() -> f64 {
      let mismatch_fraction: f64 = 2.0 / 20.0;
      -0.75 * (1.0 - 4.0 / 3.0 * mismatch_fraction).ln()
    }

    fn manual_indel_count_on_edge(reconstruction: &PartitionSet, edge_key: GraphEdgeKey) -> usize {
      sparse_reconstruction(reconstruction)
        .partition
        .obs_edges
        .get(edge_key)
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
      k.to_f64_exact() * lambda.ln() - lambda - ln_factorial(u64::try_from(k).unwrap())
    }

    pub(super) fn manual_total_indel_log_lh(
      graph: &Graph,
      reconstruction: &PartitionSet,
      branch_lengths: &EdgeMap<Option<f64>>,
    ) -> f64 {
      let total_indels: usize = graph
        .get_edges()
        .map(|edge_ref| manual_indel_count_on_edge(reconstruction, edge_ref.key()))
        .sum();
      let total_branch_length: f64 = graph
        .get_edges()
        .map(|edge_ref| branch_lengths.get(edge_ref.key()).copied().flatten().unwrap_or(0.0))
        .sum();
      let indel_rate = if total_indels > 0 && total_branch_length > 0.0 {
        total_indels.to_f64_exact() / total_branch_length
      } else {
        0.0
      };

      graph
        .get_edges()
        .map(|edge_ref| {
          let edge_key = edge_ref.key();
          let branch_length = branch_lengths.get(edge_key).copied().flatten().unwrap_or(0.0);
          let indel_count = manual_indel_count_on_edge(reconstruction, edge_key);
          manual_poisson_indel_log_lh(indel_count, indel_rate, branch_length)
        })
        .sum()
    }
  }
}
