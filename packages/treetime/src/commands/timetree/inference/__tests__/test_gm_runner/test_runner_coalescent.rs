#[cfg(test)]
mod tests {
  use super::super::test_gm_runner_support::support::{
    ALPHABET, OUTPUTS, load_alignment_for_dataset, load_dates_for_dataset,
  };
  use crate::commands::ancestral::marginal::initialize_marginal;
  use crate::commands::clock::clock_model::ClockModel;
  use crate::commands::clock::clock_regression::{ClockParams, estimate_clock_model_with_reroot};
  use crate::commands::clock::date_constraints::load_date_constraints;
  use crate::commands::clock::find_best_root::params::BranchPointOptimizationParams;
  use crate::commands::timetree::inference::runner::run_timetree;
  use crate::commands::timetree::partition_ops::PartitionTimetreeAll;
  use crate::commands::timetree::utils::{
    extract_node_times, initialize_clock_totals_from_time_distributions, initialize_node_divergences,
  };
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::representation::partition::marginal_dense::PartitionMarginalDense;
  use crate::representation::partition::timetree::GraphTimetree;
  use crate::representation::payload::timetree::{EdgeTimetree, NodeTimetree};
  use eyre::Report;
  use maplit::btreemap;
  use parking_lot::RwLock;
  use rstest::rstest;
  use std::sync::Arc;
  use treetime_distribution::Distribution;
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
  #[ignore]
  fn test_runner_coalescent_completes(#[case] tc: f64) -> Result<(), Report> {
    let dataset = "flu_h3n2_20";
    let case = &OUTPUTS[dataset];

    let (graph, partitions, clock_model) = build_timetree_setup(dataset, case)?;
    let mut graph = graph;
    let coalescent_tc = Distribution::constant(tc);
    run_timetree(&mut graph, &partitions, &clock_model, Some(&coalescent_tc))?;

    let times = extract_node_times(&graph);
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

  fn build_timetree_setup(
    dataset: &str,
    case: &super::super::test_gm_runner_support::support::DatasetOutputs,
  ) -> Result<
    (
      GraphTimetree,
      Vec<Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>>>,
      ClockModel,
    ),
    Report,
  > {
    let mut graph: GraphTimetree = nwk_read_str(case.rerooted_tree_nwk())?;
    let dates = load_dates_for_dataset(dataset)?;
    load_date_constraints(&dates, &graph)?;

    let aln = load_alignment_for_dataset(dataset)?;
    let dense_partition = Arc::new(RwLock::new(PartitionMarginalDense {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: ALPHABET.clone(),
      length: case.sequence_length(),
      nodes: btreemap! {},
      edges: btreemap! {},
    }));

    let partitions: Vec<Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>>> = vec![dense_partition];
    initialize_marginal(&graph, &partitions, &aln)?;
    initialize_node_divergences(&graph);

    let clock_model = estimate_clock_model_with_reroot(
      &mut graph,
      &ClockParams::default(),
      Some(case.clock_rate()),
      true,
      &BranchPointOptimizationParams::default(),
      None,
    )?;

    initialize_clock_totals_from_time_distributions(&graph)?;

    Ok((graph, partitions, clock_model))
  }
}
