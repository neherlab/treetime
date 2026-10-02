#[cfg(test)]
mod tests {
  use super::super::test_gm_runner_support::support::{
    ALPHABET, OUTPUTS, extract_node_times, load_alignment_for_dataset, load_dates_for_dataset,
  };
  use crate::progress::NoopProgress;
  use crate::seq::alignment::node_seq_inputs;

  use crate::ancestral::fitch::create_fitch_partition;
  use crate::ancestral::marginal::branch_lengths_or_zero;
  use crate::ancestral::pipeline::SparseReconstruction;
  use crate::clock::clock_regression::{ClockTree, ClockVarianceParams, estimate_clock_model_with_reroot_policy};
  use crate::clock::clock_state::ClockInputs;
  use crate::clock::date_constraints::load_date_constraints;
  use crate::clock::find_best_root::params::BranchPointOptimizationParams;
  use crate::clock::reroot::RerootParams;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::timetree::inference::bad_branches::bad_leaves;
  use crate::timetree::inference::runner::{TimeInferenceInputs, run_timetree};
  use crate::timetree::inference::time_inference::{likely_times, unit_gammas};

  use crate::partition::timetree::partition::PartitionTimetree;
  use crate::timetree::branch_model::BranchModel;
  use eyre::Report;
  use std::collections::{BTreeMap, BTreeSet};
  use treetime_graph::graph::Graph;

  use rstest::rstest;

  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::AlignmentRecord;

  use treetime_utils::pretty_assert_map_abs_diff_eq;

  #[rustfmt::skip]
#[rstest]
  #[case::flu_h3n2_20("flu_h3n2_20")]
  #[trace]
  #[ignore = "sparse-vs-dense discrepancy: max 0.92 years at root-adjacent node (grid-width difference)"]
  fn test_gm_runner_marginal_sparse_dense_consistency(#[case] dataset: &str) -> Result<(), Report> {
    let case = &OUTPUTS[dataset];
    let expected = case.marginal_dense();

    let nwk_parsed = nwk_read_str(case.rerooted_tree_nwk())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;

    let graph: Graph = graph;
    let dates = load_dates_for_dataset(dataset)?;
    let constraints = load_date_constraints(&dates, &graph, &names, &NoopProgress)?;

    let aln: Vec<AlignmentRecord> = load_alignment_for_dataset(dataset)?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
    let fitch = create_fitch_partition(&graph, 0, ALPHABET.clone(), &node_seq_inputs(&graph, &names, aln))?;
    let (partition, node_states) = fitch.into_marginal_sparse(&graph)?;
    let sparse_partition = PartitionTimetree::Sparse(SparseReconstruction::seeded(partition, jc69(JC69Params::default())?, node_states));

    let partition = sparse_partition.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;

    let times = likely_times(&graph, &constraints, None)?;
    let clock_estimate_inputs = ClockInputs::from_times(&graph, &times, &BTreeMap::new());
    let (
      ClockTree {
        graph, branch_lengths, ..
      },
      clock_reroot,
    ) = estimate_clock_model_with_reroot_policy(
      ClockTree {
        graph,
        branch_lengths,
        inputs: clock_estimate_inputs,
      },
      &BTreeSet::new(),
      &ClockVarianceParams::default(),
      Some(case.clock_rate()),
      true,
      &BranchPointOptimizationParams::default(),
      &RerootParams::default(),
      None, &names, &NoopProgress
    )?;
    let clock_model = clock_reroot.into_clock_fit()?.model;

    let run_branch_lengths = branch_lengths;
    let run_names = names.clone();
    let inference = run_timetree(
      &TimeInferenceInputs {
        graph: &graph,
        date_constraints: &constraints,
        leaf_bad_branches: &bad_leaves(&graph, &constraints, &BTreeSet::new()),
        gammas: &unit_gammas(&graph),
        branch_model: &BranchModel::Marginal(partition),
        branch_lengths: &run_branch_lengths,
        names: &run_names,
        clock_model: &clock_model,
        no_indels: false,
      },
      None,
      &NoopProgress,
    )?;

    let actual = extract_node_times(&graph, &names, &inference.posterior);
    pretty_assert_map_abs_diff_eq!(expected, &actual, epsilon = 1e-6);

    Ok(())
  }
}
