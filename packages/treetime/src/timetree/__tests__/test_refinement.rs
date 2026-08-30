#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::marginal::initialize_marginal;
  use crate::clock::clock_model::ClockModel;
  use crate::clock::clock_regression::{ClockParams, estimate_clock_model_with_reroot};
  use crate::clock::date_constraints::load_date_constraints;
  use crate::clock::find_best_root::params::BranchPointOptimizationParams;
  use crate::coalescent::coalescent::CoalescentModel;
  use crate::coalescent::lineage_counts::compute_lineage_counts;
  use crate::coalescent::total_lh::compute_coalescent_total_lh;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::timetree::{GraphTimetree, PartitionTimetree, PartitionTimetreeAllVec};
  use crate::payload::clock_set::ClockSet;
  use crate::pretty_assert_abs_diff_eq;
  use crate::seq::alignment::get_common_length;
  use crate::timetree::inference::runner::run_timetree;
  use crate::timetree::refinement::{
    Refinement, RefinementOptions, RefinementOutcome, TopologyOutcome, TopologyRefinement,
  };
  use crate::timetree::utils::{initialize_clock_totals_from_time_distributions, initialize_node_divergences};
  use eyre::Report;
  use indoc::indoc;
  use maplit::btreemap;
  use ndarray::{Array1, array};
  use parking_lot::RwLock;
  use pretty_assertions::assert_eq;
  use std::sync::Arc;
  use treetime_distribution::Distribution;
  use treetime_graph::edge::{BranchDistribution, HasBranchLength};
  use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
  use treetime_io::dates_csv::{DateConstraint, DatesMap};
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::assert_error;
  use treetime_utils::io::json::{JsonPretty, json_write_str};
  use treetime_utils::sync::random::get_random_number_generator;

  const CLOCK_RATE: f64 = 0.001;

  #[test]
  #[ignore = "mass-sized node times break downstream invariants (positional log-lh, polytomy resolution): kb/issues/H-timetree-mass-sizing-node-times-break-downstream-invariants.md"]
  fn test_refinement_rebuilds_complete_coalescent_state_after_topology_change() -> Result<(), Report> {
    let (mut graph, partitions, mut clock_model) = create_polytomy_state()?;
    let tc = Distribution::constant(10.0);

    let outcome = refine(&mut graph, &partitions, &mut clock_model, Some(&tc))?;

    assert_eq!(0, outcome.sequence_changes);
    assert_eq!(TopologyOutcome::Changed { resolved_nodes: 1 }, outcome.topology);
    assert_eq!(5, graph.get_nodes().len());
    assert!(
      graph
        .get_nodes()
        .iter()
        .all(|node| { node.read_arc().payload().read_arc().time_distribution.is_some() })
    );

    let edge_lh = compute_coalescent_total_lh(&graph, &tc)?;
    let model = CoalescentModel::new(&compute_lineage_counts(&graph)?, &tc)?;
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
    pretty_assert_abs_diff_eq!(node_lh, edge_lh.value(), epsilon = 1e-10);

    let outcome = refine(&mut graph, &partitions, &mut clock_model, Some(&tc))?;
    assert_eq!(TopologyOutcome::Unchanged, outcome.topology);

    Ok(())
  }

  #[test]
  fn test_refinement_missing_time_preserves_inference_state() -> Result<(), Report> {
    let (mut graph, partitions, mut clock_model) = create_polytomy_state()?;
    let root = graph.get_exactly_one_root()?;
    root.write_arc().payload().write_arc().time = None;
    let expected_error = format!(
      "Polytomy resolution failed: Polytomy resolution requires an inferred time for node {}, but it has none",
      root.read_arc().key()
    );
    let before = serialize_state(&graph, &partitions, &clock_model)?;

    assert_error!(
      refine(
        &mut graph,
        &partitions,
        &mut clock_model,
        Some(&Distribution::constant(10.0)),
      ),
      expected_error
    );

    let after = serialize_state(&graph, &partitions, &clock_model)?;
    assert_eq!(before, after);

    Ok(())
  }

  #[test]
  fn test_refinement_non_finite_time_preserves_inference_state() -> Result<(), Report> {
    let (mut graph, partitions, mut clock_model) = create_polytomy_state()?;
    let root = graph.get_exactly_one_root()?;
    let root_key = root.read_arc().key();
    root.write_arc().payload().write_arc().time = Some(f64::NAN);
    let before = serialize_state(&graph, &partitions, &clock_model)?;

    assert_error!(
      refine(
        &mut graph,
        &partitions,
        &mut clock_model,
        Some(&Distribution::constant(10.0)),
      ),
      format!(
        "Polytomy resolution failed: Polytomy resolution requires a finite inferred time for node {root_key}, but it has NaN"
      )
    );

    let after = serialize_state(&graph, &partitions, &clock_model)?;
    assert_eq!(before, after);
    assert_eq!(
      f64::NAN.to_bits(),
      root
        .read_arc()
        .payload()
        .read_arc()
        .time
        .expect("Root time must remain present")
        .to_bits()
    );

    Ok(())
  }

  #[test]
  #[ignore = "mass-sized node times break downstream invariants (positional log-lh, polytomy resolution): kb/issues/H-timetree-mass-sizing-node-times-break-downstream-invariants.md"]
  fn test_refinement_unchanged_topology_recomputes_missing_time() -> Result<(), Report> {
    let (mut graph, partitions, mut clock_model) = create_polytomy_state()?;
    let tc = Distribution::constant(10.0);
    refine(&mut graph, &partitions, &mut clock_model, Some(&tc))?;
    let root = graph.get_exactly_one_root()?;
    root.write_arc().payload().write_arc().time = None;

    let outcome = refine(&mut graph, &partitions, &mut clock_model, None)?;

    assert_eq!(TopologyOutcome::Unchanged, outcome.topology);
    assert!(root.read_arc().payload().read_arc().time.is_some_and(f64::is_finite));

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

  /// Polytomy resolution samples; pin the stream so refinement tests stay deterministic.
  const REFINEMENT_TEST_SEED: u64 = 0xC0FFEE;

  /// Timescale behind the merger rate polytomy resolution is driven with. The pipeline
  /// estimates one from the tree when a run carries no coalescent prior; these tests pin it, so
  /// the rate is the same however `coalescent_tc` is set.
  const REFINEMENT_TEST_TC: f64 = 10.0;

  fn refine(
    graph: &mut GraphTimetree,
    partitions: &PartitionTimetreeAllVec,
    clock_model: &mut ClockModel,
    coalescent_tc: Option<&Distribution>,
  ) -> Result<RefinementOutcome, Report> {
    let pinned_tc = Distribution::constant(REFINEMENT_TEST_TC);
    let coalescent = CoalescentModel::new(&compute_lineage_counts(graph)?, coalescent_tc.unwrap_or(&pinned_tc))?;
    let merger_rate =
      coalescent.branch_merger_rate_schedule(&PiecewiseConstantFn::new(array![], array![REFINEMENT_TEST_TC]))?;

    Refinement {
      graph,
      partitions,
      clock_model,
      clock_params: &ClockParams::default(),
      branch_params: &BranchPointOptimizationParams::default(),
      merger_rate: &merger_rate,
      prior: coalescent_tc.is_some().then_some(&coalescent),
      rng: &mut get_random_number_generator(Some(REFINEMENT_TEST_SEED)),
      options: &refinement_options(),
    }
    .run()
  }

  fn refinement_options() -> RefinementOptions {
    RefinementOptions {
      relax: vec![],
      topology: TopologyRefinement::Resolve,
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
