#[cfg(test)]
mod tests {
  use super::super::test_gm_runner_support::support::{
    ALPHABET, OUTPUTS, load_alignment_for_dataset, load_dates_for_dataset,
  };
  use crate::ancestral::marginal::profile_branch_lengths;
  use crate::ancestral::pipeline::DenseReconstruction;
  use crate::clock::clock_regression::{ClockParams, estimate_clock_model_with_reroot_policy};
  use crate::clock::clock_state::{ClockInputs, ClockState};
  use crate::clock::date_constraints::load_date_constraints;
  use crate::clock::find_best_root::params::BranchPointOptimizationParams;
  use crate::clock::reroot::RerootParams;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::timetree::marginal::initialize_marginal_timetree;
  use crate::partition::timetree::partition::PartitionTimetree;
  use crate::timetree::inference::runner::run_timetree;
  use crate::timetree::timetree_state::TimetreeState;
  use crate::timetree::utils::{extract_node_times, initialize_node_divergences};
  use eyre::Report;
  use treetime_graph::graph::Graph;

  use rstest::rstest;
  use treetime_io::nwk::{NwkParse, nwk_read_str};
  use treetime_utils::pretty_assert_map_abs_diff_eq;

  // --- Marginal dense tests ---

  #[rustfmt::skip]
#[rstest]
  // #[case::dengue_20("dengue_20")]       // TODO: missing internal node times, leaf dates not refined
  // #[case::ebola_20("ebola_20")]         // TODO: golden master node key mismatch (v0 captures 11 internal nodes, v1 rerooting produces 19)
  #[case::flu_h3n2_20("flu_h3n2_20")]
  // #[case::lassa_l_20("lassa_L_20")]     // TODO: missing internal node times, leaf dates not refined
  // #[case::mpox_clade_ii_20("mpox_clade_ii_20")] // TODO: missing internal node times, leaf dates not refined
  // #[case::rsv_a_20("rsv_a_20")]         // TODO: missing internal node times, leaf dates not refined
  // #[case::tb_20("tb_20")]               // TODO: missing internal node times, leaf dates not refined
  // #[case::zika_20("zika_20")]           // TODO: read_dates strips # from headers, name_column="#name" mismatches
  #[trace]
  #[ignore = "dense-vs-v0 discrepancy: max 0.92 years at root-adjacent node (grid-width difference)"]
  // TODO: investigate why uniform branch distribution grids produce 0.92-year shift at
  // root-adjacent nodes. All other nodes within 0.1 years. Related:
  // kb/issues/M-timetree-branch-grid-uniform-resolution.md
  fn test_gm_runner_marginal_dense(#[case] dataset: &str) -> Result<(), Report> {
    let case = &OUTPUTS[dataset];
    let expected = case.marginal_dense();

    let NwkParse { graph, names, mut branch_lengths, .. } = nwk_read_str(case.rerooted_tree_nwk())?;

    let mut graph: Graph = graph;
    let dates = load_dates_for_dataset(dataset)?;
    let constraints = load_date_constraints(&dates, &graph, &names)?;

    let aln = load_alignment_for_dataset(dataset)?;
    let dense_partition = PartitionTimetree::Dense(DenseReconstruction {
      partition: PartitionMarginalDense::new(0, jc69(JC69Params::default())?, ALPHABET.clone(), case.sequence_length()),
      node_states: std::collections::BTreeMap::new(),
      backward: std::collections::BTreeMap::new(),
      forward: std::collections::BTreeMap::new(),
      estimates: std::collections::BTreeMap::new(),
    });

    let mut partitions: Vec<PartitionTimetree> = vec![dense_partition];
    initialize_marginal_timetree(&graph, &profile_branch_lengths(&branch_lengths), &mut partitions, &aln, &names)?.value();
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
      None, &names_tt_1
    )?;
    let clock_model = clock_reroot.into_clock_model()?;

    let mut state = TimetreeState::new(&graph);
    let run_branch_lengths = branch_lengths;
    let run_names = names.clone();
    state = run_timetree(
      &mut graph,
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
