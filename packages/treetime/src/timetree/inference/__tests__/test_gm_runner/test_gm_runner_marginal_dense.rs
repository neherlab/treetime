#[cfg(test)]
mod tests {
  use super::super::test_gm_runner_support::support::{
    ALPHABET, OUTPUTS, extract_node_times, load_alignment_for_dataset, load_dates_for_dataset,
  };
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::clock::clock_regression::{ClockTree, ClockVarianceParams, estimate_clock_model_with_reroot_policy};
  use crate::clock::clock_state::ClockInputs;
  use crate::clock::date_constraints::load_date_constraints;
  use crate::clock::find_best_root::params::BranchPointOptimizationParams;
  use crate::clock::reroot::RerootParams;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::marginal::reconstruction::DenseReconstruction;
  use crate::partition::marginal::reconstruction::MarginalReconstruction;
  use crate::progress::NoopProgress;
  use crate::seq::alignment::node_seq_inputs;
  use crate::test_utils::dates_by_node;
  use crate::timetree::branch_model::BranchModel;
  use crate::timetree::inference::bad_branches::bad_leaves;
  use crate::timetree::inference::result::given_times;
  use crate::timetree::inference::runner::{TimeInferenceInputs, run_timetree};
  use crate::timetree::optimization::relaxed_clock::unit_gammas;
  use eyre::Report;
  use rstest::rstest;
  use std::collections::{BTreeMap, BTreeSet};
  use treetime_grid::MaxGridPoints;
  use treetime_io::nwk::nwk_read;
  use treetime_primitives::AlignmentRecord;
  use treetime_utils::pretty_assert_map_abs_diff_eq;

  #[rustfmt::skip]
  #[rstest]
  #[case::flu_h3n2_20("flu_h3n2_20")]
  // TODO: enable these datasets when their golden-master gaps are fixed: kb/issues/M-timetree-gm-runner-missing-internal-times.md,
  // kb/issues/M-timetree-dense-golden-master-node-mismatch.md, kb/issues/M-timetree-date-header-hash.md
  // #[case::dengue_20("dengue_20")]       // TODO: missing internal node times, leaf dates not refined
  // #[case::ebola_20("ebola_20")]         // TODO: golden master node key mismatch (v0 captures 11 internal nodes, v1 rerooting produces 19)
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

    let nwk_parsed = nwk_read(case.rerooted_tree_nwk().as_bytes())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;

    let dates = load_dates_for_dataset(dataset)?;
    let constraints = load_date_constraints(&dates_by_node(dates, &graph, &names), &graph, &NoopProgress)?;

    let aln: Vec<AlignmentRecord> = load_alignment_for_dataset(dataset)?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
    let dense_partition = MarginalReconstruction::Dense(DenseReconstruction::seeded(
      PartitionMarginalDense::new(0, ALPHABET.clone(), &graph, &node_seq_inputs(&graph, &names, aln))?,
      jc69(JC69Params::default())?,
    ));

    let (partition, _) = dense_partition.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;

    let times = given_times(&graph, &constraints)?;
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

    let run_names = names.clone();
    let inference = run_timetree(
      &TimeInferenceInputs {
        graph: &graph,
        date_constraints: &constraints,
        leaf_bad_branches: &bad_leaves(&graph, &constraints, &BTreeSet::new()),
        gammas: &unit_gammas(&graph),
        branch_model: &BranchModel::Marginal(partition),
        branch_lengths: &branch_lengths,
        names: &run_names,
        clock_model: &clock_model,
        clock_rate_fixed: false,
        no_indels: false,
        max_grid_points: MaxGridPoints::default(),
      },
      None,
      &NoopProgress,
    )?;

    let actual = extract_node_times(&graph, &names, &inference.posterior);
    pretty_assert_map_abs_diff_eq!(expected, &actual, epsilon = 1e-6);

    Ok(())
  }
}
