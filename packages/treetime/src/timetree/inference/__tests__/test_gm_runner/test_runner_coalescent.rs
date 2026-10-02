#[cfg(test)]
mod tests {
  use super::super::test_gm_runner_support::support::{
    ALPHABET, OUTPUTS, extract_node_times, load_alignment_for_dataset, load_dates_for_dataset,
  };
  use crate::ancestral::marginal::branch_lengths_or_zero;
  use crate::clock::clock_model::ClockModel;
  use crate::clock::clock_regression::{ClockTree, ClockVarianceParams, estimate_clock_model_with_reroot_policy};
  use crate::clock::clock_state::ClockInputs;
  use crate::clock::date_constraints::{DateConstraints, load_date_constraints};
  use crate::clock::find_best_root::params::BranchPointOptimizationParams;
  use crate::clock::reroot::RerootParams;
  use crate::coalescent::coalescent::CoalescentModel;
  use crate::coalescent::lineage_counts::compute_lineage_counts;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::marginal::reconstruction::DenseReconstruction;
  use crate::partition::marginal::reconstruction::MarginalReconstruction;
  use crate::partition::marginal::shared::update::MarginalEdges;
  use crate::progress::NoopProgress;
  use crate::seq::alignment::node_seq_inputs;
  use crate::test_utils::constraint_coalescent_node_times;
  use crate::timetree::branch_model::BranchModel;
  use crate::timetree::inference::bad_branches::bad_leaves;
  use crate::timetree::inference::runner::{TimeInferenceInputs, run_timetree};
  use crate::timetree::inference::time_inference::{likely_times, unit_gammas};
  use eyre::Report;
  use helpers::build_timetree_setup;
  use treetime_graph::graph::Graph;

  use rstest::rstest;
  use std::collections::{BTreeMap, BTreeSet};
  use treetime_distribution::Distribution;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::AlignmentRecord;

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

let (graph, names, partition, clock_model, constraints, branch_lengths) = build_timetree_setup(dataset, case)?;
    let node_times = constraint_coalescent_node_times(&graph, &constraints)?;
    let coalescent = CoalescentModel::new(&compute_lineage_counts(&graph, &node_times)?, &Distribution::constant(tc))?;
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
      Some(&coalescent),
      &NoopProgress,
    )?;

    let times = extract_node_times(&graph, &names, &inference.posterior);
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
    MarginalReconstruction,
    ClockModel,
    DateConstraints,
    BTreeMap<GraphEdgeKey, Option<f64>>,
  );

  mod helpers {
    use super::*;

    pub(super) fn build_timetree_setup(
      dataset: &str,
      case: &super::super::super::test_gm_runner_support::support::DatasetOutputs,
    ) -> Result<TimetreeSetup, Report> {
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
      let dense_partition = MarginalReconstruction::Dense(DenseReconstruction {
        partition: PartitionMarginalDense::new(0, ALPHABET.clone(), &graph, &node_seq_inputs(&graph, &names, aln))?,
        gtr: jc69(JC69Params::default())?,
        node_states: BTreeMap::new(),
        edges: MarginalEdges::default(),
      });

      let (partition, _) = dense_partition.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;

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
        None,
        &names,
        &NoopProgress,
      )?;
      let clock_model = clock_reroot.into_clock_fit()?.model;

      Ok((graph, names, partition, clock_model, constraints, branch_lengths))
    }
  }
}
