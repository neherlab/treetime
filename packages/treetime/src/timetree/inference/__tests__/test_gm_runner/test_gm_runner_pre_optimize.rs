#[cfg(test)]
mod tests {
  use super::super::test_gm_runner_support::support::{
    ALPHABET, OUTPUTS, load_alignment_for_dataset, load_dates_for_dataset,
  };
  use crate::seq::alignment::node_seq_inputs;

  use crate::ancestral::fitch::create_fitch_partition;
  use crate::ancestral::marginal::branch_lengths_or_zero;
  use crate::ancestral::pipeline::SparseReconstruction;
  use crate::clock::clock_regression::{ClockVarianceParams, estimate_clock_model_with_reroot_policy};
  use crate::clock::clock_state::{ClockInputs, ClockState};
  use crate::clock::date_constraints::load_date_constraints;
  use crate::clock::find_best_root::params::BranchPointOptimizationParams;
  use crate::clock::reroot::RerootParams;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::dispatch::run_optimize_mixed;
  use crate::optimize::gather::{
    gather_timetree_edge_contributions, gather_timetree_edge_indel_counts, timetree_total_sequence_length,
  };
  use crate::optimize::params::BranchOptMethod;
  use crate::partition::timetree::marginal::{initialize_marginal_timetree, marginal_update_timetree};
  use crate::timetree::inference::runner::run_timetree;
  use crate::timetree::timetree_state::TimetreeState;
  use crate::timetree::utils::{extract_node_times, initialize_node_divergences};
  use std::collections::BTreeMap;
  use treetime_graph::edge::GraphEdgeKey;

  use crate::partition::timetree::partition::PartitionTimetree;
  use eyre::Report;
  use itertools::Itertools;
  use treetime_graph::graph::Graph;

  use rstest::rstest;

  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::AlignmentRecord;

  fn extract_branch_lengths(graph: &Graph, branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>) -> Vec<f64> {
    graph
      .get_edges()
      .map(|e| branch_lengths[&e.key()].unwrap_or(0.0))
      .collect_vec()
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::flu_h3n2_20("flu_h3n2_20")]
  #[trace]
  fn test_gm_runner_pre_optimize_changes_branch_lengths(#[case] dataset: &str) -> Result<(), Report> {
    let case = &OUTPUTS[dataset];

    let nwk_parsed = nwk_read_str(case.rerooted_tree_nwk())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;

    let graph: Graph = graph;
    let aln: Vec<AlignmentRecord> = load_alignment_for_dataset(dataset)?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
    let fitch = create_fitch_partition(&graph, 0, ALPHABET.clone(), &node_seq_inputs(&graph, &names, aln.clone()))?;
    let (partition, node_states) = fitch.into_marginal_sparse(&graph)?;
    let sparse_partition = PartitionTimetree::Sparse(SparseReconstruction::seeded(partition, jc69(JC69Params::default())?, node_states));

    let partitions: Vec<PartitionTimetree> = vec![sparse_partition];
    let (partitions, _) = initialize_marginal_timetree(&graph, &branch_lengths_or_zero(&branch_lengths), partitions, &node_seq_inputs(&graph, &names, aln))?;

    let before = extract_branch_lengths(&graph, &branch_lengths);

    let total_length = timetree_total_sequence_length(&partitions);
    let contributions = gather_timetree_edge_contributions(&graph, &partitions)?;
    let indel_counts = gather_timetree_edge_indel_counts(&graph, &partitions);
    run_optimize_mixed(
      &graph,
      total_length,
      &contributions,
      &indel_counts,
      BranchOptMethod::BrentSqrt,
      &mut branch_lengths,
    )?;

    let after = extract_branch_lengths(&graph, &branch_lengths);

    let n_changed = before
      .iter()
      .zip(&after)
      .filter(|(b, a)| (*b - *a).abs() > 1e-10)
      .count();

    assert!(
      n_changed > 0,
      "Expected at least one branch length to change after ML optimization, but none did"
    );

    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::flu_h3n2_20("flu_h3n2_20")]
  #[trace]
  fn test_gm_runner_pre_optimize_pipeline_succeeds(#[case] dataset: &str) -> Result<(), Report> {
    let case = &OUTPUTS[dataset];

    let nwk_parsed = nwk_read_str(case.rerooted_tree_nwk())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;

    let mut graph: Graph = graph;
    let dates = load_dates_for_dataset(dataset)?;
    let constraints = load_date_constraints(&dates, &graph, &names)?;

    let aln: Vec<AlignmentRecord> = load_alignment_for_dataset(dataset)?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
    let fitch = create_fitch_partition(&graph, 0, ALPHABET.clone(), &node_seq_inputs(&graph, &names, aln.clone()))?;
    let (partition, node_states) = fitch.into_marginal_sparse(&graph)?;
    let sparse_partition = PartitionTimetree::Sparse(SparseReconstruction::seeded(partition, jc69(JC69Params::default())?, node_states));

    let partitions: Vec<PartitionTimetree> = vec![sparse_partition];
    let (partitions, _) = initialize_marginal_timetree(&graph, &branch_lengths_or_zero(&branch_lengths), partitions, &node_seq_inputs(&graph, &names, aln))?;
    let mut clock_state = ClockState::new(&graph);
    initialize_node_divergences(&graph, &mut clock_state, &branch_lengths, &names)?;

    let total_length = timetree_total_sequence_length(&partitions);
    let contributions = gather_timetree_edge_contributions(&graph, &partitions)?;
    let indel_counts = gather_timetree_edge_indel_counts(&graph, &partitions);
    run_optimize_mixed(
      &graph,
      total_length,
      &contributions,
      &indel_counts,
      BranchOptMethod::BrentSqrt,
      &mut branch_lengths,
    )?;
    let (partitions, _) = marginal_update_timetree(&graph, &branch_lengths_or_zero(&branch_lengths), partitions)?;

    let times = TimetreeState::seed_from_values(&graph, &constraints).likely_times(&constraints);
    let mut clock_estimate_inputs = ClockInputs::seed_from_times(&graph, &times);
    let names_tt_1 = names.clone();
    let clock_estimate_state = ClockState::new(&graph);
    let (_clock_estimate_state, clock_reroot) = estimate_clock_model_with_reroot_policy(
      &mut graph,
      &mut clock_estimate_inputs,
      clock_estimate_state,
      &ClockVarianceParams::default(),
      Some(case.clock_rate()),
      true,
      &BranchPointOptimizationParams::default(),
      &RerootParams::default(),
      &mut branch_lengths,
      None, &names_tt_1
    )?;
    let clock_model = clock_reroot.into_clock_model()?;

    let mut state = TimetreeState::seed_from_values(&graph, &constraints);
    let run_branch_lengths = branch_lengths;
    let run_names = names.clone();
    state = run_timetree(
      &graph,
      &constraints,
      &partitions,      &run_branch_lengths,
      &run_names,
      &clock_model,
      None,
      false,
      state,
      &mut clock_state,
    )?;

    let actual = extract_node_times(&graph, &names, &state);
    assert!(
      !actual.is_empty(),
      "Expected node times to be populated after timetree inference with pre-optimization"
    );

    for (name, time) in &actual {
      assert!(
        time.is_finite(),
        "Node {name} has non-finite time {time} after timetree inference with pre-optimization"
      );
    }

    Ok(())
  }
}
