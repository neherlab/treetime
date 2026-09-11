#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::marginal::initialize_marginal;
  use crate::ancestral::marginal::profile_branch_lengths;
  use crate::clock::clock_model::ClockModel;
  use crate::clock::clock_regression::{ClockParams, estimate_clock_model_with_reroot_policy};
  use crate::clock::clock_state::ClockState;
  use crate::clock::date_constraints::load_date_constraints;
  use crate::clock::find_best_root::params::BranchPointOptimizationParams;
  use crate::clock::reroot::RerootParams;
  use crate::coalescent::coalescent::CoalescentModel;
  use crate::coalescent::lineage_counts::compute_lineage_counts;
  use crate::coalescent::total_lh::compute_coalescent_total_lh;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::timetree::partition::{GraphTimetree, PartitionTimetree};
  use crate::pretty_assert_abs_diff_eq;
  use crate::seq::alignment::get_common_length;
  use crate::timetree::inference::runner::run_timetree;
  use crate::timetree::refinement::{
    Refinement, RefinementOptions, RefinementOutcome, TopologyOutcome, TopologyRefinement,
  };
  use crate::timetree::timetree_state::TimetreeState;
  use crate::timetree::utils::initialize_node_divergences;
  use eyre::Report;
  use indoc::indoc;
  use maplit::btreemap;
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime_distribution::Distribution;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::node::GraphNodeKey;
  use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
  use treetime_io::dates_csv::{DateConstraint, DatesMap};
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::{NwkParse, nwk_read_str};
  use treetime_utils::assert_error;
  use treetime_utils::io::json::{JsonPretty, json_write_str};
  use treetime_utils::sync::random::get_random_number_generator;

  const CLOCK_RATE: f64 = 0.001;

  #[test]
  #[ignore = "mass-sized node times break downstream invariants (positional log-lh, polytomy resolution): kb/issues/H-timetree-mass-sizing-node-times-break-downstream-invariants.md"]
  fn test_refinement_rebuilds_complete_coalescent_state_after_topology_change() -> Result<(), Report> {
    let (mut graph, names, mut partitions, mut clock_model, mut state, mut branch_lengths) = create_polytomy_state()?;
    let tc = Distribution::constant(10.0);

    let outcome = refine(
      &mut graph,
      &names,
      &mut partitions,
      &mut clock_model,
      Some(&tc),
      &mut state,
      &mut branch_lengths,
    )?;

    assert_eq!(0, outcome.sequence_changes);
    assert_eq!(TopologyOutcome::Changed { resolved_nodes: 1 }, outcome.topology);
    assert_eq!(5, graph.get_nodes().len());
    assert!(
      graph
        .get_nodes()
        .iter()
        .all(|node| { state.node(node.read_arc().key()).time_distribution.is_some() })
    );

    let edge_lh = compute_coalescent_total_lh(&graph, &tc, &state.coalescent_node_times())?;
    let model = CoalescentModel::new(&compute_lineage_counts(&graph, &state.coalescent_node_times())?, &tc)?;
    let node_lh = -graph
      .get_nodes()
      .iter()
      .map(|node| {
        let node = node.read_arc();
        let time = state
          .node(node.key())
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

    let outcome = refine(
      &mut graph,
      &names,
      &mut partitions,
      &mut clock_model,
      Some(&tc),
      &mut state,
      &mut branch_lengths,
    )?;
    assert_eq!(TopologyOutcome::Unchanged, outcome.topology);

    Ok(())
  }

  #[test]
  fn test_refinement_missing_time_preserves_inference_state() -> Result<(), Report> {
    let (mut graph, names, mut partitions, mut clock_model, mut state, mut branch_lengths) = create_polytomy_state()?;
    let root_key = graph.get_exactly_one_root()?.read_arc().key();
    state.node_mut(root_key).time = None;
    let expected_error = format!(
      "Polytomy resolution failed: Polytomy resolution requires an inferred time for node {root_key}, but it has none"
    );
    let before = serialize_state(&graph, &partitions, &clock_model)?;

    assert_error!(
      refine(
        &mut graph,
        &names,
        &mut partitions,
        &mut clock_model,
        Some(&Distribution::constant(10.0)),
        &mut state,
        &mut branch_lengths,
      ),
      expected_error
    );

    let after = serialize_state(&graph, &partitions, &clock_model)?;
    assert_eq!(before, after);

    Ok(())
  }

  #[test]
  fn test_refinement_non_finite_time_preserves_inference_state() -> Result<(), Report> {
    let (mut graph, names, mut partitions, mut clock_model, mut state, mut branch_lengths) = create_polytomy_state()?;
    let root_key = graph.get_exactly_one_root()?.read_arc().key();
    state.node_mut(root_key).time = Some(f64::NAN);
    let before = serialize_state(&graph, &partitions, &clock_model)?;

    assert_error!(
      refine(
        &mut graph,
        &names,
        &mut partitions,
        &mut clock_model,
        Some(&Distribution::constant(10.0)),
        &mut state,
        &mut branch_lengths,
      ),
      format!(
        "Polytomy resolution failed: Polytomy resolution requires a finite inferred time for node {root_key}, but it has NaN"
      )
    );

    let after = serialize_state(&graph, &partitions, &clock_model)?;
    assert_eq!(before, after);
    assert_eq!(
      f64::NAN.to_bits(),
      state
        .node(root_key)
        .time
        .expect("Root time must remain present")
        .to_bits()
    );

    Ok(())
  }

  #[test]
  #[ignore = "mass-sized node times break downstream invariants (positional log-lh, polytomy resolution): kb/issues/H-timetree-mass-sizing-node-times-break-downstream-invariants.md"]
  fn test_refinement_unchanged_topology_recomputes_missing_time() -> Result<(), Report> {
    let (mut graph, names, mut partitions, mut clock_model, mut state, mut branch_lengths) = create_polytomy_state()?;
    let tc = Distribution::constant(10.0);
    refine(
      &mut graph,
      &names,
      &mut partitions,
      &mut clock_model,
      Some(&tc),
      &mut state,
      &mut branch_lengths,
    )?;
    let root_key = graph.get_exactly_one_root()?.read_arc().key();
    state.node_mut(root_key).time = None;

    let outcome = refine(
      &mut graph,
      &names,
      &mut partitions,
      &mut clock_model,
      None,
      &mut state,
      &mut branch_lengths,
    )?;

    assert_eq!(TopologyOutcome::Unchanged, outcome.topology);
    assert!(state.node(root_key).time.is_some_and(f64::is_finite));

    Ok(())
  }

  fn create_polytomy_state() -> Result<
    (
      GraphTimetree,
      BTreeMap<GraphNodeKey, Option<String>>,
      Vec<PartitionTimetree>,
      ClockModel,
      TimetreeState,
      BTreeMap<GraphEdgeKey, Option<f64>>,
    ),
    Report,
  > {
    let NwkParse {
      graph,
      names,
      mut branch_lengths,
      ..
    } = nwk_read_str("(A:0.01,B:0.01,C:0.01)root;")?;
    let mut graph: GraphTimetree = graph;
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
    let mut partitions = vec![PartitionTimetree::Dense(PartitionMarginalDense::new(
      0,
      jc69(JC69Params::default())?,
      alphabet,
      get_common_length(&aln)?,
    ))];
    initialize_marginal(
      &graph,
      &profile_branch_lengths(&branch_lengths),
      &mut partitions,
      &aln,
      &names,
    )?;

    let dates: DatesMap = btreemap! {
      "A".to_owned() => Some(DateConstraint::exact(2010.0)),
      "B".to_owned() => Some(DateConstraint::exact(2015.0)),
      "C".to_owned() => Some(DateConstraint::exact(2020.0)),
    };
    let constraints = load_date_constraints(&dates, &graph, &names)?;
    let mut clock_state = ClockState::new(&graph);
    initialize_node_divergences(&graph, &mut clock_state, &branch_lengths, &names)?;

    let times = TimetreeState::seed_from_values(&graph, &constraints).likely_times();
    let mut clock_estimate_state = ClockState::seed_from_values(&graph, &times);
    let names_tt_1 = names.clone();
    let clock_model = estimate_clock_model_with_reroot_policy(
      &mut graph,
      &mut clock_estimate_state,
      &ClockParams::default(),
      Some(CLOCK_RATE),
      true,
      &BranchPointOptimizationParams::default(),
      &RerootParams::default(),
      &mut branch_lengths,
      None,
      &names_tt_1,
    )?
    .into_clock_model()?;
    let mut state = TimetreeState::seed_from_values(&graph, &constraints);
    let run_branch_lengths = branch_lengths;
    let run_names = names.clone();
    run_timetree(
      &mut graph,
      &mut partitions,
      &run_branch_lengths,
      &run_names,
      &clock_model,
      None,
      false,
      &mut state,
      &mut clock_state,
    )?;

    Ok((graph, names, partitions, clock_model, state, run_branch_lengths))
  }

  fn serialize_state(
    graph: &GraphTimetree,
    partitions: &[PartitionTimetree],
    clock_model: &ClockModel,
  ) -> Result<SerializedState, Report> {
    let graph = json_write_str(graph, JsonPretty(false))?;
    let partitions = partitions
      .iter()
      .map(|partition| json_write_str(&partition, JsonPretty(false)))
      .collect::<Result<Vec<_>, _>>()?;
    let clock_model = json_write_str(clock_model, JsonPretty(false))?;
    Ok(SerializedState {
      graph,
      partitions,
      clock_model,
    })
  }

  /// Polytomy resolution samples; pin the stream so refinement tests stay deterministic.
  const REFINEMENT_TEST_SEED: u64 = 0xC0FFEE;

  /// Timescale behind the merger rate polytomy resolution is driven with. The pipeline
  /// estimates one from the tree when a run carries no coalescent prior; these tests pin it, so
  /// the rate is the same however `coalescent_tc` is set.
  const REFINEMENT_TEST_TC: f64 = 10.0;

  fn refine(
    graph: &mut GraphTimetree,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    partitions: &mut [PartitionTimetree],
    clock_model: &mut ClockModel,
    coalescent_tc: Option<&Distribution>,
    state: &mut TimetreeState,
    branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  ) -> Result<RefinementOutcome, Report> {
    let pinned_tc = Distribution::constant(REFINEMENT_TEST_TC);
    let coalescent = CoalescentModel::new(
      &compute_lineage_counts(graph, &state.coalescent_node_times())?,
      coalescent_tc.unwrap_or(&pinned_tc),
    )?;
    let merger_rate =
      coalescent.branch_merger_rate_schedule(&PiecewiseConstantFn::new(array![], array![REFINEMENT_TEST_TC]))?;

    let mut clock_state = ClockState::new(graph);
    let mut clock_branch_lengths: BTreeMap<GraphEdgeKey, f64> = BTreeMap::new();
    let mut names = names.clone();

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
      state,
      clock_state: &mut clock_state,
      clock_branch_lengths: &mut clock_branch_lengths,
      branch_lengths,
      names: &mut names,
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
  }
}
