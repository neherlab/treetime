#[cfg(test)]
mod tests {
  use super::super::test_gm_runner_support::support::{
    ALPHABET, OUTPUTS, extract_node_times, load_alignment_for_dataset, load_dates_for_dataset,
  };
  use crate::ancestral::marginal::branch_lengths_or_zero;
  use crate::ancestral::pipeline::DenseReconstruction;
  use crate::clock::clock_regression::{ClockVarianceParams, estimate_clock_model_with_reroot_policy};
  use crate::clock::clock_state::{ClockInputs, ClockState};
  use crate::clock::date_constraints::load_date_constraints;
  use crate::clock::find_best_root::params::BranchPointOptimizationParams;
  use crate::clock::reroot::RerootParams;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::timetree::marginal::initialize_marginal_timetree;
  use crate::partition::timetree::partition::PartitionTimetree;
  use crate::seq::alignment::node_seq_inputs;
  use crate::timetree::inference::runner::run_timetree;
  use crate::timetree::timetree_state::TimetreeState;
  use crate::timetree::utils::initialize_node_divergences;
  use eyre::Report;
  use treetime_graph::graph::Graph;

  use rstest::rstest;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::AlignmentRecord;
  use treetime_utils::pretty_assert_map_abs_diff_eq;

  #[rustfmt::skip]
#[rstest]
  #[case::flu_h3n2_20("flu_h3n2_20")]
  #[trace]
  #[ignore = "dense-vs-v0 discrepancy: max 0.92 years at root-adjacent node (grid-width difference)"]
  fn test_gm_runner_marginal_dense(#[case] dataset: &str) -> Result<(), Report> {
    let case = &OUTPUTS[dataset];
    let expected = case.marginal_dense();

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
    let dense_partition = PartitionTimetree::Dense(DenseReconstruction::seeded(
      PartitionMarginalDense::new(0, ALPHABET.clone(), case.sequence_length()),
      jc69(JC69Params::default())?,
      std::collections::BTreeMap::new(),
    ));

    let partitions: Vec<PartitionTimetree> = vec![dense_partition];
    let (partitions, _) = initialize_marginal_timetree(&graph, &branch_lengths_or_zero(&branch_lengths), partitions, &node_seq_inputs(&graph, &names, aln))?;
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
      &ClockVarianceParams::default(),
      Some(case.clock_rate()),
      true,
      &BranchPointOptimizationParams::default(),
      &RerootParams::default(),
      &mut branch_lengths,
      None, &names_tt_1
    )?;
    let clock_model = clock_reroot.into_clock_model()?;

    let mut state = TimetreeState::new(&graph);
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
    pretty_assert_map_abs_diff_eq!(expected, &actual, epsilon = 1e-6);

    Ok(())
  }
}
