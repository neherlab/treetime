#[cfg(test)]
mod tests {
  use treetime_io::nwk::nwk_fasta_node_inputs;
  use super::super::test_gm_runner_support::support::{
    ALPHABET, OUTPUTS, load_alignment_for_dataset, load_dates_for_dataset,
  };
  use crate::ancestral::marginal::profile_branch_lengths;
  use crate::ancestral::pipeline::DenseReconstruction;
  use crate::clock::clock_model::ClockModel;
  use crate::clock::clock_regression::{ClockParams, estimate_clock_model_with_reroot_policy};
  use crate::clock::clock_state::{ClockInputs, ClockState};
  use crate::clock::date_constraints::{DateConstraints, load_date_constraints};
  use crate::clock::find_best_root::params::BranchPointOptimizationParams;
  use crate::clock::reroot::RerootParams;
  use crate::coalescent::coalescent::CoalescentModel;
  use crate::coalescent::lineage_counts::compute_lineage_counts;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::marginal::shared::update::MarginalEdges;
  use crate::partition::timetree::marginal::initialize_marginal_timetree;
  use crate::partition::timetree::partition::PartitionTimetree;
  use crate::timetree::inference::runner::run_timetree;
  use crate::timetree::timetree_state::TimetreeState;
  use crate::timetree::utils::{extract_node_times, initialize_node_divergences};
  use eyre::Report;
  use treetime_graph::graph::Graph;

  use rstest::rstest;
  use std::collections::BTreeMap;
  use treetime_distribution::Distribution;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::nwk::nwk_read_str;

  /// Verifies the full timetree pipeline completes without panic when coalescent
  /// is enabled. Before the Formula discretization fixes, this would panic on
  /// Formula * Function multiplication or on Formula methods like likely_time().
  #[rustfmt::skip]
  #[rstest]
  #[case::tc_0_1( 0.1)]
  #[case::tc_1_0( 1.0)]
  #[case::tc_10_0(10.0)]
  #[trace]
  #[ignore = "golden master datasets not yet passing"]
  fn test_runner_coalescent_completes(#[case] tc: f64) -> Result<(), Report> {
    let dataset = "flu_h3n2_20";
    let case = &OUTPUTS[dataset];

let (graph, names, partitions, clock_model, constraints, branch_lengths) = build_timetree_setup(dataset, case)?;
    let mut graph = graph;
    let node_times = TimetreeState::seed_from_values(&graph, &constraints).coalescent_node_times();
    let coalescent = CoalescentModel::new(&compute_lineage_counts(&graph, &node_times)?, &Distribution::constant(tc))?;
    let mut state = TimetreeState::new(&graph);
    let mut clock_state = ClockState::new(&graph);
    let run_branch_lengths = branch_lengths;
    let run_names = names.clone();
    state = run_timetree(
      &mut graph,
      &constraints,
      &partitions,      &run_branch_lengths,
      &run_names,
      &clock_model,
      Some(&coalescent),
      false,
      state,
      &mut clock_state,
    )?;

    let times = extract_node_times(&graph, &names, &state);
    let expected_count = graph.num_nodes();
    assert_eq!(
      expected_count,
      times.len(),
      "All nodes (tips + internal) must have time values, but only {actual} of {expected} do",
      actual = times.len(),
      expected = expected_count,
    );

    for (name, time) in &times {
      assert!(time.is_finite(), "Node {name} has non-finite time {time}");
    }

    Ok(())
  }

  type TimetreeSetup = (
    Graph,
    BTreeMap<GraphNodeKey, Option<String>>,
    Vec<PartitionTimetree>,
    ClockModel,
    DateConstraints,
    BTreeMap<GraphEdgeKey, Option<f64>>,
  );

  fn build_timetree_setup(
    dataset: &str,
    case: &super::super::test_gm_runner_support::support::DatasetOutputs,
  ) -> Result<TimetreeSetup, Report> {
    let nwk_parsed = nwk_read_str(case.rerooted_tree_nwk())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let dates = load_dates_for_dataset(dataset)?;
    let constraints = load_date_constraints(&dates, &graph, &names)?;

    let aln = load_alignment_for_dataset(dataset)?;
    let dense_partition = PartitionTimetree::Dense(DenseReconstruction {
      partition: PartitionMarginalDense::new(
        0,
        jc69(JC69Params::default())?,
        ALPHABET.clone(),
        case.sequence_length(),
      ),
      node_states: BTreeMap::new(),
      edges: MarginalEdges::default(),
    });

    let partitions: Vec<PartitionTimetree> = vec![dense_partition];
    let (partitions, _) = initialize_marginal_timetree(&graph, &profile_branch_lengths(&branch_lengths), partitions, &nwk_fasta_node_inputs(&graph, &names, aln))?;
    let mut clock_state = ClockState::new(&graph);
    initialize_node_divergences(&graph, &mut clock_state, &branch_lengths, &names)?;

    let times = TimetreeState::seed_from_values(&graph, &constraints).likely_times(&constraints);
    let mut clock_estimate_inputs = ClockInputs::seed_from_times(&graph, &times);
    let names_tt_1 = names.clone();
    let clock_estimate_state = ClockState::new(&graph);
    let (_clock_estimate_state, clock_reroot) = estimate_clock_model_with_reroot_policy(
      &mut graph,
      &mut clock_estimate_inputs,
      clock_estimate_state,
      &ClockParams::default(),
      Some(case.clock_rate()),
      true,
      &BranchPointOptimizationParams::default(),
      &RerootParams::default(),
      &mut branch_lengths,
      None,
      &names_tt_1,
    )?;
    let clock_model = clock_reroot.into_clock_model()?;

    Ok((graph, names, partitions, clock_model, constraints, branch_lengths))
  }
}
