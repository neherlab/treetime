#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::branch_lengths::{branch_lengths_or_zero, one_mutation};
  use crate::clock::clock_regression::{ClockFit, ClockFitParams, ClockTree, fit_clock_model};
  use crate::clock::clock_state::ClockInputs;
  use crate::clock::date_constraints::{DateConstraints, load_date_constraints};
  use crate::coalescent::mean_field::{CoalescentPrior, LineageExpectation};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::marginal::reconstruction::DenseReconstruction;
  use crate::partition::marginal::shared::update::MarginalEdges;
  use crate::partition::set::PartitionSet;
  use crate::pretty_assert_abs_diff_eq;
  use crate::pretty_assert_ulps_eq;
  use crate::progress::NoopProgress;
  use crate::test_utils::dates_by_node;
  use crate::test_utils::jc69_input_branch_model;
  use crate::test_utils::leaf_seq_inputs;
  use crate::test_utils::{find_node_key_by_name, marginal_timetree_params, parent_edge_key};
  use crate::timetree::branch_model::BranchModel;
  use crate::timetree::convergence::sequence_changes::capture_ancestral_states;
  use crate::timetree::inference::bad_branches::bad_leaves;
  use crate::timetree::inference::result::{PosteriorSummaries, given_times};
  use crate::timetree::inference::runner::{TimeInferenceInputs, run_timetree};
  use crate::timetree::optimization::marginal_clock_rate::{MarginalRateInputs, RateWindow, step_marginal_clock_rate};
  use crate::timetree::optimization::relaxed_clock::{apply_relaxed_clock, unit_gammas};
  use crate::timetree::params::{BranchLengthMode, FinalSequences, TimetreeContext, TimetreeParams};
  use crate::timetree::round::{RoundInputs, RoundOutcome, RoundState, TopologyOutcome, refinement_round};
  use eyre::Report;
  use indoc::indoc;
  use itertools::Itertools;
  use maplit::btreemap;
  use treetime_graph::key_map::EdgeMap;

  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::{BTreeMap, BTreeSet};
  use treetime_distribution::Distribution;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::node::GraphNodeKey;
  use treetime_graph::schedule::GraphSchedule;
  use treetime_grid::MaxGridPoints;

  use treetime_io::dates_csv::DateConstraint;
  use treetime_io::fasta::fasta_read;
  use treetime_io::nwk::{NwkNodeComments, NwkWriteOptions, nwk_read, nwk_write_str};
  use treetime_primitives::AlignmentRecord;
  use treetime_utils::assert_error;
  use treetime_utils::sync::random::get_random_number_generator;

  const CLOCK_RATE: f64 = 0.001;

  const ROUND_TEST_SEED: u64 = 0xC0FFEE;

  const ROUND_TEST_TC: f64 = 10.0;

  const UNARY_BRANCH_LENGTH: f64 = 0.01;

  const UNARY_TREE: &str = "((A:0.01,B:0.01,C:0.01)P:0.01,(D:0.01)S:0.01)root;";

  const UNARY_FASTA: &str = indoc! {r#"
    >A
    ACGTACGTACGTACGTACGTACGT
    >B
    ACGTACGAACGTACGTACGTACGT
    >C
    ACGTACGTACGTACCTACGTACGT
    >D
    ACGTTCGTACGTACGTACGTACGA
  "#};

  const TIME_AFTER_EVERY_SAMPLE: f64 = 2100.0;

  const RELAX: [f64; 2] = [1.0, 1.0];

  const RATE_TREE: &str = "((A:0.018,B:0.024,C:0.030)F:0.012,(D:0.006,E:0.010)S:0.004)root;";

  const RATE_TEST_CLOCK_RATE: f64 = 0.002;

  const RATE_SEQUENCE_LENGTH: usize = 1000;

  #[test]
  fn test_round_rebuilds_complete_coalescent_state_after_topology_change() -> Result<(), Report> {
    let (state, context) = helpers::create_polytomy_state()?;
    let tc = Distribution::constant(10.0);

    let (state, outcome) = helpers::refine(&context, state, Some(&tc))?;

    assert_eq!(0, outcome.sequence_changes);
    assert_eq!(1, outcome.topology.resolved_nodes());
    assert_eq!(5, state.graph.get_nodes().count());
    let inference = &state.time_inference;
    assert!(
      state
        .graph
        .get_nodes()
        .all(|node| inference.posterior[&node.key()].likely_time.is_some())
    );

    let (_, outcome) = helpers::refine(&context, state, Some(&tc))?;
    assert_eq!(0, outcome.topology.resolved_nodes());

    Ok(())
  }

  #[test]
  fn test_round_missing_time_fails_polytomy_resolution() -> Result<(), Report> {
    let (mut state, context) = helpers::create_polytomy_state()?;
    let root_key = state.graph.get_exactly_one_root()?.key();
    state
      .time_inference
      .posterior
      .get_mut(root_key)
      .expect("root must have a posterior")
      .time = None;

    assert_error!(
      helpers::refine(&context, state, Some(&Distribution::constant(10.0))),
      format!(
        "Polytomy resolution failed: Polytomy resolution requires an inferred time for node {root_key}, but it has none"
      )
    );

    Ok(())
  }

  #[test]
  fn test_round_non_finite_time_fails_polytomy_resolution() -> Result<(), Report> {
    let (mut state, context) = helpers::create_polytomy_state()?;
    let root_key = state.graph.get_exactly_one_root()?.key();
    state
      .time_inference
      .posterior
      .get_mut(root_key)
      .expect("root must have a posterior")
      .time = Some(f64::NAN);

    assert_error!(
      helpers::refine(&context, state, Some(&Distribution::constant(10.0))),
      format!(
        "Polytomy resolution failed: Polytomy resolution requires a finite inferred time for node {root_key}, but it has NaN"
      )
    );

    Ok(())
  }

  #[test]
  fn test_round_unchanged_topology_recomputes_missing_time() -> Result<(), Report> {
    let (state, context) = helpers::create_polytomy_state()?;
    let tc = Distribution::constant(10.0);
    let (mut state, _) = helpers::refine(&context, state, Some(&tc))?;
    let root_key = state.graph.get_exactly_one_root()?.key();
    state
      .time_inference
      .posterior
      .get_mut(root_key)
      .expect("root must have a posterior")
      .time = None;

    let (state, outcome) = helpers::refine(&context, state, None)?;

    assert_eq!(0, outcome.topology.resolved_nodes());
    assert!(
      state.time_inference.posterior[&root_key]
        .time
        .is_some_and(f64::is_finite)
    );

    Ok(())
  }

  #[test]
  fn test_round_is_deterministic_for_equal_inputs_and_seeds() -> Result<(), Report> {
    let (state_a, context) = helpers::create_polytomy_state()?;
    let (state_b, _) = helpers::create_polytomy_state()?;

    let (state_a, outcome_a) = helpers::refine(&context, state_a, None)?;
    let (state_b, outcome_b) = helpers::refine(&context, state_b, None)?;

    assert_eq!(state_a.time_inference, state_b.time_inference);
    assert_eq!(state_a.names, state_b.names);
    assert_eq!(state_a.branch_lengths, state_b.branch_lengths);
    assert_eq!(state_a.clock_branch_lengths, state_b.clock_branch_lengths);
    assert_eq!(state_a.gammas, state_b.gammas);
    assert_eq!(state_a.clock_points, state_b.clock_points);
    assert_eq!(helpers::newick(&state_a)?, helpers::newick(&state_b)?);
    assert_eq!(
      state_a.clock_model.clock_rate().to_bits(),
      state_b.clock_model.clock_rate().to_bits()
    );
    assert_eq!(
      state_a.clock_model.intercept().to_bits(),
      state_b.clock_model.intercept().to_bits()
    );
    assert_eq!(helpers::outcome_parts(&outcome_a), helpers::outcome_parts(&outcome_b));
    Ok(())
  }

  #[test]
  fn test_round_polytomy_round_that_only_removes_a_single_child_node_reports_a_topology_change() -> Result<(), Report> {
    let (state, context) = helpers::create_unary_state()?;
    let state = helpers::without_room_above_the_polytomy(state);

    let (state, outcome) = helpers::refine(&context, state, None)?;

    assert!(matches!(
      &outcome.topology,
      TopologyOutcome::Changed { merger_times } if merger_times.is_empty()
    ));
    let node_names: BTreeSet<String> = state
      .graph
      .get_nodes()
      .map(|node| state.names[node.key()].to_owned())
      .collect();
    let expected_names: BTreeSet<String> = ["A", "B", "C", "D", "P", "root"].map(str::to_owned).into();
    assert_eq!(expected_names, node_names);
    helpers::assert_state_matches_graph(&state);
    let d_key = find_node_key_by_name(&state.graph, &state.names, "D").expect("leaf D must remain");
    let merged = state.branch_lengths[&parent_edge_key(&state.graph, d_key)].expect("the merged edge keeps a length");
    pretty_assert_abs_diff_eq!(merged, 2.0 * UNARY_BRANCH_LENGTH, epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_round_single_child_removal_gives_the_merged_edge_the_fitted_gamma_of_the_child_edge() -> Result<(), Report> {
    let (state, context) = helpers::create_unary_state()?;
    let state = helpers::without_room_above_the_polytomy(state);
    let fitted = helpers::fitted_gammas(&context, &state)?;
    let d_key = find_node_key_by_name(&state.graph, &state.names, "D").expect("leaf D must exist");
    let child_edge_gamma = fitted[&parent_edge_key(&state.graph, d_key)];

    let (state, outcome) = helpers::refine_with(&context, state, None, &helpers::round_params(true, RELAX.to_vec()))?;

    assert!(outcome.topology.resolved_nodes() == 0 && outcome.topology.changed());
    let merged_edge = parent_edge_key(&state.graph, d_key);
    let expected: EdgeMap<f64> = state
      .graph
      .get_edges()
      .map(|edge| {
        let gamma = if edge.key() == merged_edge {
          child_edge_gamma
        } else {
          fitted[&edge.key()]
        };
        (edge.key(), gamma)
      })
      .collect();
    assert_eq!(expected, state.gammas);
    assert!(
      expected.values().any(|gamma| (gamma - 1.0).abs() > 1e-3),
      "the relaxed clock must move at least one kept rate away from 1: {expected:?}"
    );
    Ok(())
  }

  #[test]
  fn test_round_polytomy_resolution_keeps_fitted_gammas_and_starts_new_edges_at_one() -> Result<(), Report> {
    let (state, context) = helpers::create_wide_polytomy_state()?;
    let fitted = helpers::fitted_gammas(&context, &state)?;

    let (state, outcome) = helpers::refine_with(&context, state, None, &helpers::round_params(true, RELAX.to_vec()))?;

    let TopologyOutcome::Changed { merger_times } = &outcome.topology else {
      panic!("the round must resolve the polytomy");
    };
    assert!(!merger_times.is_empty());
    let expected: EdgeMap<f64> = state
      .graph
      .get_edges()
      .map(|edge| {
        let gamma = if merger_times.contains_key(edge.target()) {
          1.0
        } else {
          fitted[&edge.key()]
        };
        (edge.key(), gamma)
      })
      .collect();
    assert_eq!(expected, state.gammas);
    assert!(
      expected.values().any(|gamma| (gamma - 1.0).abs() > 1e-3),
      "the relaxed clock must move at least one kept rate away from 1: {expected:?}"
    );
    Ok(())
  }

  #[test]
  fn test_round_polytomy_resolution_keeps_a_dated_single_child_node_at_its_date() -> Result<(), Report> {
    let (state, context) = helpers::create_state(
      UNARY_TREE,
      UNARY_FASTA,
      &btreemap! {
        "A".to_owned() => Some(DateConstraint::exact(2010.0)),
        "B".to_owned() => Some(DateConstraint::exact(2012.0)),
        "C".to_owned() => Some(DateConstraint::exact(2014.0)),
        "D".to_owned() => Some(DateConstraint::exact(2016.0)),
        "S".to_owned() => Some(DateConstraint::exact(2008.0)),
      },
    )?;
    let state = helpers::without_room_above_the_polytomy(state);

    let (state, _) = helpers::refine(&context, state, None)?;

    let s_key = find_node_key_by_name(&state.graph, &state.names, "S").expect("dated node S must remain");
    assert_eq!(
      1,
      state.graph.get_node(s_key).expect("Node must exist").outbound().len()
    );
    helpers::assert_state_matches_graph(&state);
    let time = state.time_inference.posterior[&s_key].time.expect("S must have a time");
    pretty_assert_abs_diff_eq!(2008.0, time, epsilon = 1e-9);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::input(   BranchLengthMode::Input)]
  #[case::marginal(BranchLengthMode::Marginal)]
  #[trace]
  fn test_round_relaxed_clock_fits_rates_to_the_durations_between_node_dates(
    #[case] mode: BranchLengthMode,
  ) -> Result<(), Report> {
    let (state, context) = helpers::create_rate_state(mode)?;
    let params = TimetreeParams {
      clock_rate: Some(RATE_TEST_CLOCK_RATE),
      ..helpers::round_params(true, RELAX.to_vec())
    };
    let first_fit = helpers::fitted_gammas(&context, &state)?;

    let (state, outcome) = helpers::refine_with(&context, state, None, &params)?;

    let TopologyOutcome::Changed { merger_times } = &outcome.topology else {
      panic!("the first round must resolve the polytomy");
    };
    let expected: EdgeMap<f64> = state
      .graph
      .get_edges()
      .map(|edge| {
        let gamma = if merger_times.contains_key(edge.target()) {
          1.0
        } else {
          first_fit[&edge.key()]
        };
        (edge.key(), gamma)
      })
      .collect();
    assert_eq!(expected, state.gammas);

    let second_fit = helpers::fitted_gammas(&context, &state)?;
    let (state, outcome) = helpers::refine_with(&context, state, None, &params)?;

    assert!(!outcome.topology.changed());
    assert_eq!(second_fit, state.gammas);
    let leaf_gamma = |name: &str| {
      let key = find_node_key_by_name(&state.graph, &state.names, name).expect("leaf must exist");
      state.gammas[&parent_edge_key(&state.graph, key)]
    };
    let slowest_fast = ["A", "B", "C"].map(leaf_gamma).into_iter().fold(f64::INFINITY, f64::min);
    let fastest_slow = ["D", "E"].map(leaf_gamma).into_iter().fold(f64::NEG_INFINITY, f64::max);
    assert!(
      slowest_fast > fastest_slow,
      "every edge of the clade that evolves 3 times faster must get a higher rate than every edge of the slow clade: \
       slowest fast {slowest_fast}, fastest slow {fastest_slow}"
    );
    Ok(())
  }

  #[test]
  fn test_round_unchanged_topology_keeps_the_relaxed_gammas() -> Result<(), Report> {
    let (state, context) = helpers::create_unary_state()?;
    let expected = helpers::fitted_gammas(&context, &state)?;

    let (state, outcome) = helpers::refine_with(&context, state, None, &helpers::round_params(false, RELAX.to_vec()))?;

    assert!(!outcome.topology.changed());
    assert_eq!(expected, state.gammas);
    assert!(
      state.gammas.values().any(|gamma| (gamma - 1.0).abs() > 1e-3),
      "the relaxed clock must move at least one rate away from 1: {:?}",
      state.gammas
    );
    Ok(())
  }

  #[test]
  fn test_round_rate_step_after_a_topology_change_uses_the_rebuilt_coalescent_prior() -> Result<(), Report> {
    let (state, context) = helpers::create_state(
      "((A:0.001,B:0.009,C:0.009,D:0.017)P:0.004,E:0.017)root;",
      indoc! {r#"
        >A
        ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
        >B
        ACGTACGAACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
        >C
        ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACCTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
        >D
        ACGTACGTACGTACGTACGTTCGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
        >E
        TCGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAAGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
      "#},
      &btreemap! {
        "A".to_owned() => Some(DateConstraint::exact(2010.0)),
        "B".to_owned() => Some(DateConstraint::exact(2012.0)),
        "C".to_owned() => Some(DateConstraint::exact(2014.0)),
        "D".to_owned() => Some(DateConstraint::exact(2016.0)),
        "E".to_owned() => Some(DateConstraint::exact(2015.0)),
      },
    )?;
    let reference_rate = state.clock_model.clock_rate();
    let expectation = LineageExpectation::of_inference(&state.graph, &state.time_inference)?;
    let prior = CoalescentPrior::new(&expectation, &Distribution::constant(ROUND_TEST_TC))?;
    let params = TimetreeParams {
      clock_rate: None,
      ..helpers::round_params(true, vec![])
    };

    let (state, outcome) = helpers::refine_with_prior(&context, state, Some(&prior), &params)?;

    let TopologyOutcome::Changed { merger_times } = &outcome.topology else {
      panic!("the round must resolve the polytomy");
    };
    assert!(!merger_times.is_empty());
    let rate_step = |prior: &CoalescentPrior| {
      step_marginal_clock_rate(&MarginalRateInputs {
        graph: &state.graph,
        date_constraints: &context.date_constraints,
        coalescent: Some(prior),
        bad_branches: &state.time_inference.bad_branches,
        branches: &state.time_inference.branches,
        reference_rate,
        window: RateWindow::initial(),
        max_points: MaxGridPoints::default(),
      })
      .map(|step| step.clock_rate)
    };
    let rebuilt_prior_rate = rate_step(&prior.with_topology(&state.graph, merger_times)?)?;
    assert_ne!(rebuilt_prior_rate.to_bits(), rate_step(&prior)?.to_bits());
    pretty_assert_ulps_eq!(rebuilt_prior_rate, state.clock_model.clock_rate(), max_ulps = 0);
    Ok(())
  }

  mod helpers {
    use super::*;
    use crate::optimize::indel::IndelPolicy;
    use pretty_assertions::assert_eq;
    use treetime_graph::key_map::EdgeMap;
    use treetime_graph::key_map::NodeMap;

    pub(super) fn newick(state: &RoundState) -> Result<String, Report> {
      nwk_write_str(
        &GraphSchedule::new(&state.graph)?.into_tree()?,
        &state.names,
        &state.branch_lengths,
        &NwkWriteOptions::default(),
        &NwkNodeComments::new(),
      )
    }

    pub(super) fn outcome_parts(outcome: &RoundOutcome) -> (usize, usize, Option<u64>, Option<u64>) {
      (
        outcome.sequence_changes,
        outcome.topology.resolved_nodes(),
        outcome.time_change.max.map(f64::to_bits),
        outcome.time_change.rms.map(f64::to_bits),
      )
    }

    pub(super) fn create_polytomy_state() -> Result<(RoundState, TimetreeContext), Report> {
      create_state(
        "(A:0.01,B:0.01,C:0.01)root;",
        indoc! {r#"
          >A
          ACGTACGTACGT
          >B
          ACGTACGTACGT
          >C
          ACGTACGTACGT
        "#},
        &btreemap! {
          "A".to_owned() => Some(DateConstraint::exact(2010.0)),
          "B".to_owned() => Some(DateConstraint::exact(2015.0)),
          "C".to_owned() => Some(DateConstraint::exact(2020.0)),
        },
      )
    }

    pub(super) fn create_unary_state() -> Result<(RoundState, TimetreeContext), Report> {
      create_state(
        UNARY_TREE,
        UNARY_FASTA,
        &btreemap! {
          "A".to_owned() => Some(DateConstraint::exact(2010.0)),
          "B".to_owned() => Some(DateConstraint::exact(2012.0)),
          "C".to_owned() => Some(DateConstraint::exact(2014.0)),
          "D".to_owned() => Some(DateConstraint::exact(2016.0)),
        },
      )
    }

    pub(super) fn create_wide_polytomy_state() -> Result<(RoundState, TimetreeContext), Report> {
      create_state(
        "((A:0.01,B:0.08,C:0.06,D:0.005)P:0.02,(E:0.03,F:0.01)Y:0.01)root;",
        indoc! {r#"
          >A
          ACGTACGTACGTACGTACGTACGT
          >B
          ACGTACGAACGTACGTACGTACGT
          >C
          ACGTACGTACGTACCTACGTACGT
          >D
          ACGTTCGTACGTACGTACGTACGA
          >E
          TCGTACGTACGTACGTACGTACGT
          >F
          ACGTACGTACGAACGTACGTACGT
        "#},
        &btreemap! {
          "A".to_owned() => Some(DateConstraint::exact(2010.0)),
          "B".to_owned() => Some(DateConstraint::exact(2012.0)),
          "C".to_owned() => Some(DateConstraint::exact(2014.0)),
          "D".to_owned() => Some(DateConstraint::exact(2016.0)),
          "E".to_owned() => Some(DateConstraint::exact(2011.0)),
          "F".to_owned() => Some(DateConstraint::exact(2015.0)),
        },
      )
    }

    pub(super) fn fitted_gammas(context: &TimetreeContext, state: &RoundState) -> Result<EdgeMap<f64>, Report> {
      apply_relaxed_clock(
        &state.graph,
        &state.branch_lengths,
        &RELAX,
        one_mutation(context.sequence_length),
        state.clock_model.clock_rate(),
        &node_date_durations(state),
      )
    }

    fn node_date_durations(state: &RoundState) -> EdgeMap<Option<f64>> {
      let times = state.time_inference.node_times();
      state
        .graph
        .get_edges()
        .map(|edge| {
          let duration = times[&edge.source()]
            .zip(times[&edge.target()])
            .map(|(parent, child)| child - parent);
          (edge.key(), duration)
        })
        .collect()
    }

    pub(super) fn create_rate_state(mode: BranchLengthMode) -> Result<(RoundState, TimetreeContext), Report> {
      create_state_with(
        RATE_TREE,
        &rate_fasta(),
        &btreemap! {
          "A".to_owned() => Some(DateConstraint::exact(2010.0)),
          "B".to_owned() => Some(DateConstraint::exact(2012.0)),
          "C".to_owned() => Some(DateConstraint::exact(2014.0)),
          "D".to_owned() => Some(DateConstraint::exact(2010.0)),
          "E".to_owned() => Some(DateConstraint::exact(2014.0)),
        },
        RATE_TEST_CLOCK_RATE,
        mode,
      )
    }

    fn rate_fasta() -> String {
      let leaves = [
        ("A", [0..12, 12..30]),
        ("B", [0..12, 30..54]),
        ("C", [0..12, 54..84]),
        ("D", [84..88, 88..94]),
        ("E", [84..88, 94..104]),
      ];
      leaves
        .iter()
        .map(|(name, mutated)| {
          let sequence: String = (0..RATE_SEQUENCE_LENGTH)
            .map(|site| {
              let bases = if mutated.iter().any(|range| range.contains(&site)) {
                b"CGTA"
              } else {
                b"ACGT"
              };
              char::from(bases[site % 4])
            })
            .collect();
          [">", name, "\n", &sequence, "\n"].concat()
        })
        .collect()
    }

    pub(super) fn create_state(
      newick: &str,
      fasta: &str,
      dates: &BTreeMap<String, Option<DateConstraint>>,
    ) -> Result<(RoundState, TimetreeContext), Report> {
      create_state_with(newick, fasta, dates, CLOCK_RATE, BranchLengthMode::Marginal)
    }

    fn create_state_with(
      newick: &str,
      fasta: &str,
      dates: &BTreeMap<String, Option<DateConstraint>>,
      clock_rate: f64,
      mode: BranchLengthMode,
    ) -> Result<(RoundState, TimetreeContext), Report> {
      let nwk_parsed = nwk_read(newick.as_bytes())?;
      let names = nwk_parsed.names.clone();
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;
      let alphabet = Alphabet::new(AlphabetName::Nuc)?;
      let aln: Vec<AlignmentRecord> = fasta_read(fasta.as_bytes(), &alphabet)?
        .into_iter()
        .map(AlignmentRecord::from)
        .collect();
      let (partition, _) = PartitionSet::Dense(vec![DenseReconstruction {
        partition: PartitionMarginalDense::new(alphabet, &graph, &leaf_seq_inputs(&graph, &names, aln))?,
        gtr: jc69(JC69Params::default())?,
        node_states: NodeMap::new(),
        edges: MarginalEdges::default(),
        layout: None,
      }])
      .marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
      let sequence_length = partition.total_length();
      let branch_model = match mode {
        BranchLengthMode::Marginal => BranchModel::Marginal(partition),
        BranchLengthMode::Input => jc69_input_branch_model(sequence_length),
      };

      let constraints = load_date_constraints(&dates_by_node(dates.clone(), &graph, &names), &graph, &NoopProgress)?;

      let times = given_times(&graph, &constraints)?;
      let inputs = ClockInputs::from_times(&graph, &times, &EdgeMap::new());
      let tree = ClockTree {
        graph,
        branch_lengths,
        inputs,
      };
      let clock_fit = fit_clock_model(
        &tree,
        &ClockFitParams {
          outliers: &BTreeSet::new(),
          clock_rate: Some(clock_rate),
          prev_clock_rate: None,
        },
        &NoopProgress,
      )?;
      let ClockTree {
        graph, branch_lengths, ..
      } = tree;
      let ClockFit {
        model: clock_model,
        points: clock_points,
      } = clock_fit.into_clock_fit()?;
      let gammas = unit_gammas(&graph);
      let time_inference = run_timetree(
        &TimeInferenceInputs {
          graph: &graph,
          date_constraints: &constraints,
          leaf_bad_branches: &bad_leaves(&graph, &constraints, &BTreeSet::new()),
          gammas: &gammas,
          branch_model: &branch_model,
          branch_lengths: &branch_lengths,
          names: &names,
          clock_model: &clock_model,
          clock_rate_fixed: false,
          indel_policy: IndelPolicy::Include,
          max_grid_points: MaxGridPoints::default(),
          posterior_summaries: PosteriorSummaries {
            regions: false,
            lineage_survival: true,
          },
        },
        None,
        &NoopProgress,
      )?;

      let state = RoundState {
        graph,
        names,
        branch_model,
        branch_lengths,
        clock_model,
        clock_points,
        rate_window: RateWindow::initial(),
        clock_branch_lengths: EdgeMap::new(),
        gammas,
        time_inference,
      };
      Ok((state, context(constraints, sequence_length)))
    }

    pub(super) fn refine(
      context: &TimetreeContext,
      state: RoundState,
      coalescent_tc: Option<&Distribution>,
    ) -> Result<(RoundState, RoundOutcome), Report> {
      refine_with(context, state, coalescent_tc, &round_params(true, vec![]))
    }

    pub(super) fn round_params(resolve_polytomies: bool, relax: Vec<f64>) -> TimetreeParams {
      TimetreeParams {
        clock_rate: Some(CLOCK_RATE),
        resolve_polytomies,
        relax,
        seed: ROUND_TEST_SEED,
        ..marginal_timetree_params()
      }
    }

    pub(super) fn refine_with(
      context: &TimetreeContext,
      state: RoundState,
      coalescent_tc: Option<&Distribution>,
      params: &TimetreeParams,
    ) -> Result<(RoundState, RoundOutcome), Report> {
      let expectation = LineageExpectation::of_inference(&state.graph, &state.time_inference)?;
      let prior = coalescent_tc
        .map(|tc| CoalescentPrior::new(&expectation, tc))
        .transpose()?;
      refine_with_prior(context, state, prior.as_ref(), params)
    }

    pub(super) fn refine_with_prior(
      context: &TimetreeContext,
      state: RoundState,
      prior: Option<&CoalescentPrior>,
      params: &TimetreeParams,
    ) -> Result<(RoundState, RoundOutcome), Report> {
      let expectation = LineageExpectation::of_inference(&state.graph, &state.time_inference)?;
      let merger_rate = expectation.merger_rate_schedule(&Distribution::constant(ROUND_TEST_TC))?;
      let leaf_bad_branches = bad_leaves(&state.graph, &context.date_constraints, &BTreeSet::new());
      let outliers = BTreeSet::new();
      let inputs = RoundInputs {
        params,
        context,
        leaf_bad_branches: &leaf_bad_branches,
        outliers: &outliers,
      };

      refinement_round(
        &inputs,
        Some(&merger_rate),
        prior,
        state,
        &mut get_random_number_generator(ROUND_TEST_SEED),
        &NoopProgress,
      )
    }

    pub(super) fn without_room_above_the_polytomy(mut state: RoundState) -> RoundState {
      let polytomy = find_node_key_by_name(&state.graph, &state.names, "P").expect("polytomy P must exist");
      state
        .time_inference
        .posterior
        .get_mut(polytomy)
        .expect("P must have a posterior")
        .time = Some(TIME_AFTER_EVERY_SAMPLE);
      state
    }

    pub(super) fn assert_state_matches_graph(state: &RoundState) {
      let nodes: BTreeSet<GraphNodeKey> = state.graph.get_nodes().map(|node| node.key()).collect();
      let internal: BTreeSet<GraphNodeKey> = state
        .graph
        .get_nodes()
        .filter(|node| !node.is_leaf())
        .map(|node| node.key())
        .collect();
      let edges: BTreeSet<GraphEdgeKey> = state.graph.get_edges().map(|edge| edge.key()).collect();
      let node_names: BTreeSet<&str> = nodes.iter().map(|&key| &state.names[key]).collect();
      assert_eq!(nodes.len(), node_names.len());
      assert_eq!(nodes, state.time_inference.posterior.keys().collect::<BTreeSet<_>>());
      assert_eq!(edges, state.gammas.keys().collect::<BTreeSet<_>>());
      assert_eq!(edges, state.branch_lengths.keys().collect::<BTreeSet<_>>());
      assert_eq!(
        internal,
        capture_ancestral_states(&state.graph, &state.branch_model)
          .iter()
          .exactly_one()
          .expect("one partition")
          .keys()
          .collect::<BTreeSet<_>>()
      );
    }

    fn context(date_constraints: DateConstraints, sequence_length: usize) -> TimetreeContext {
      TimetreeContext {
        final_sequences: FinalSequences::Skip,
        date_constraints,
        sequence_length,
      }
    }
  }
}
