#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::marginal::initialize_marginal;
  use crate::clock::clock_model::ClockModel;
  use crate::clock::clock_regression::{ClockParams, estimate_clock_model_with_reroot};
  use crate::clock::date_constraints::load_date_constraints;
  use crate::clock::find_best_root::params::BranchPointOptimizationParams;
  use crate::coalescent::coalescent::compute_coalescent_model;
  use crate::coalescent::total_lh::compute_coalescent_total_lh;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::marginal_dense::PartitionMarginalDense;
  use crate::partition::timetree::{GraphTimetree, PartitionTimetree, PartitionTimetreeAllVec};
  use crate::payload::clock_set::ClockSet;
  use crate::pretty_assert_abs_diff_eq;
  use crate::seq::alignment::get_common_length;
  use crate::timetree::inference::runner::run_timetree;
  use crate::timetree::refinement::{RefinementParams, run_refinement_iteration};
  use crate::timetree::utils::{initialize_clock_totals_from_time_distributions, initialize_node_divergences};
  use eyre::Report;
  use indoc::indoc;
  use maplit::btreemap;
  use ndarray::Array1;
  use parking_lot::RwLock;
  use pretty_assertions::assert_eq;
  use std::sync::Arc;
  use treetime_distribution::Distribution;
  use treetime_graph::edge::{BranchDistribution, HasBranchLength};
  use treetime_io::dates_csv::{DateConstraint, DatesMap};
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::assert_error;
  use treetime_utils::io::json::{JsonPretty, json_write_str};

  const CLOCK_RATE: f64 = 0.001;

  #[test]
  fn test_refinement_rebuilds_complete_coalescent_state_after_topology_change() -> Result<(), Report> {
    let (mut graph, partitions, mut clock_model) = create_polytomy_state()?;
    let tc = Distribution::constant(10.0);

    let (n_diff, n_resolved) = run_refinement_iteration(
      &refinement_params(),
      &mut graph,
      &partitions,
      &mut clock_model,
      &ClockParams::default(),
      &BranchPointOptimizationParams::default(),
      Some(&tc),
    )?;

    assert_eq!(0, n_diff);
    assert_eq!(1, n_resolved);
    assert_eq!(5, graph.get_nodes().len());
    assert!(
      graph
        .get_nodes()
        .iter()
        .all(|node| { node.read_arc().payload().read_arc().time_distribution.is_some() })
    );

    let edge_lh = compute_coalescent_total_lh(&graph, &tc)?;
    let model = compute_coalescent_model(&graph, &tc)?;
    let node_lh = -graph
      .get_nodes()
      .iter()
      .map(|node| {
        let node = node.read_arc();
        let time = node
          .payload()
          .read_arc()
          .time_distribution
          .as_ref()
          .and_then(|distribution| distribution.likely_time())
          .expect("refined node must have a likely time");
        if node.is_leaf() {
          Ok(model.leaf_contribution(time))
        } else if node.is_root() {
          model.root_contribution(time, node.outbound().len())
        } else {
          model.internal_contribution(time, node.outbound().len())
        }
      })
      .sum::<Result<f64, Report>>()?;

    // Kingman's node and edge factorizations telescope to the same objective.
    pretty_assert_abs_diff_eq!(node_lh, edge_lh, epsilon = 1e-10);

    let (_, n_resolved) = run_refinement_iteration(
      &refinement_params(),
      &mut graph,
      &partitions,
      &mut clock_model,
      &ClockParams::default(),
      &BranchPointOptimizationParams::default(),
      Some(&tc),
    )?;
    assert_eq!(0, n_resolved);

    Ok(())
  }

  #[test]
  fn test_refinement_missing_time_preserves_inference_state() -> Result<(), Report> {
    let (mut graph, partitions, mut clock_model) = create_polytomy_state()?;
    let root = graph.get_exactly_one_root()?;
    root.write_arc().payload().write_arc().time = None;
    let expected_error = format!(
      "Topology rebuild requires an inferred time for every internal node, but node {:?} has none",
      root.read_arc().key()
    );
    let before = serialize_state(&graph, &partitions, &clock_model)?;

    assert_error!(
      run_refinement_iteration(
        &refinement_params(),
        &mut graph,
        &partitions,
        &mut clock_model,
        &ClockParams::default(),
        &BranchPointOptimizationParams::default(),
        Some(&Distribution::constant(10.0)),
      ),
      expected_error
    );

    let after = serialize_state(&graph, &partitions, &clock_model)?;
    assert_eq!(before, after);

    Ok(())
  }

  fn create_polytomy_state() -> Result<(GraphTimetree, PartitionTimetreeAllVec, ClockModel), Report> {
    let mut graph: GraphTimetree = nwk_read_str("(A:0.01,B:0.01,C:0.01)root;")?;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let aln = read_many_fasta_str(
      indoc! {r#"
        >A
        ACGTACGTACGT
        >B
        ACGTACGTACGT
        >C
        ACGTACGTACGT
      "#},
      &alphabet,
    )?;
    let partitions = vec![Arc::new(RwLock::new(PartitionTimetree::Dense(
      PartitionMarginalDense::new(0, jc69(JC69Params::default())?, alphabet, get_common_length(&aln)?),
    )))];
    initialize_marginal(&graph, &partitions, &aln)?;

    let dates: DatesMap = btreemap! {
      "A".to_owned() => Some(DateConstraint::exact(2010.0)),
      "B".to_owned() => Some(DateConstraint::exact(2015.0)),
      "C".to_owned() => Some(DateConstraint::exact(2020.0)),
    };
    load_date_constraints(&dates, &graph)?;
    initialize_node_divergences(&graph)?;
    initialize_clock_totals_from_time_distributions(&graph)?;

    let clock_model = estimate_clock_model_with_reroot(
      &mut graph,
      &ClockParams::default(),
      Some(CLOCK_RATE),
      true,
      &BranchPointOptimizationParams::default(),
      None,
    )?;
    run_timetree(&mut graph, &partitions, &clock_model, None, false)?;

    let times = Array1::linspace(0.0, 30.0, 301);
    let values = times.mapv(|time: f64| (-0.5 * time).exp());
    let branch_distribution = Arc::new(Distribution::function(times, values)?);
    for edge in graph.get_edges() {
      let mut payload = edge.read_arc().payload().write_arc();
      payload.set_branch_length(Some(0.0));
      payload.set_branch_length_distribution(Some(Arc::clone(&branch_distribution)));
    }

    Ok((graph, partitions, clock_model))
  }

  fn serialize_state(
    graph: &GraphTimetree,
    partitions: &PartitionTimetreeAllVec,
    clock_model: &ClockModel,
  ) -> Result<SerializedState, Report> {
    let edge_clocks = graph_edges_clock_state(graph);
    let graph = json_write_str(graph, JsonPretty(false))?;
    let partitions = partitions
      .iter()
      .map(|partition| json_write_str(&*partition.read_arc(), JsonPretty(false)))
      .collect::<Result<Vec<_>, _>>()?;
    let clock_model = json_write_str(clock_model, JsonPretty(false))?;
    Ok(SerializedState {
      graph,
      partitions,
      clock_model,
      edge_clocks,
    })
  }

  fn graph_edges_clock_state(graph: &GraphTimetree) -> Vec<EdgeClockState> {
    graph
      .get_edges()
      .iter()
      .map(|edge| {
        let edge = edge.read_arc().payload().read_arc();
        EdgeClockState {
          to_parent: edge.clock_to_parent.clone(),
          to_child: edge.clock_to_child.clone(),
          from_child: edge.clock_from_child.clone(),
        }
      })
      .collect()
  }

  fn refinement_params() -> RefinementParams {
    RefinementParams {
      relax: vec![],
      resolve_polytomies: true,
      clock_rate: Some(CLOCK_RATE),
      no_indels: false,
    }
  }

  #[derive(Debug, PartialEq)]
  struct SerializedState {
    graph: String,
    partitions: Vec<String>,
    clock_model: String,
    edge_clocks: Vec<EdgeClockState>,
  }

  #[derive(Debug, PartialEq)]
  struct EdgeClockState {
    to_parent: ClockSet,
    to_child: ClockSet,
    from_child: ClockSet,
  }
}
