#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::clock::clock_regression::{
    ClockFit, ClockTree, ClockVarianceParams, estimate_clock_model_with_reroot_policy,
  };
  use crate::clock::clock_state::ClockInputs;
  use crate::clock::date_constraints::{DateConstraints, load_date_constraints};
  use crate::clock::find_best_root::params::BranchPointOptimizationParams;
  use crate::clock::reroot::RerootParams;
  use crate::coalescent::coalescent::CoalescentModel;
  use crate::coalescent::lineage_counts::compute_lineage_counts;
  use crate::coalescent::total_lh::compute_coalescent_total_lh;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::marginal::reconstruction::DenseReconstruction;
  use crate::partition::marginal::reconstruction::MarginalReconstruction;
  use crate::partition::marginal::shared::update::MarginalEdges;
  use crate::pretty_assert_abs_diff_eq;
  use crate::progress::NoopProgress;
  use crate::seq::alignment::node_seq_inputs;
  use crate::test_utils::marginal_timetree_params;
  use crate::timetree::branch_model::BranchModel;
  use crate::timetree::inference::bad_branches::bad_leaves;
  use crate::timetree::inference::runner::{TimeInferenceInputs, run_timetree};
  use crate::timetree::inference::time_inference::{likely_times, unit_gammas};
  use crate::timetree::params::TimeMarginalMode;
  use crate::timetree::params::{TimetreeContext, TimetreeParams};
  use crate::timetree::round::{RoundInputs, RoundOutcome, RoundState, refinement_round};
  use eyre::Report;
  use indoc::indoc;
  use maplit::btreemap;
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use std::collections::{BTreeMap, BTreeSet};
  use treetime_distribution::Distribution;
  use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
  use treetime_io::dates_csv::{DateConstraint, DatesMap};
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::{NwkWriteOptions, nwk_read_str, nwk_write_str};
  use treetime_primitives::AlignmentRecord;
  use treetime_utils::assert_error;
  use treetime_utils::sync::random::get_random_number_generator;

  const CLOCK_RATE: f64 = 0.001;

  const ROUND_TEST_SEED: u64 = 0xC0FFEE;

  const ROUND_TEST_TC: f64 = 10.0;

  #[test]
  #[ignore = "mass-sized node times break downstream invariants (positional log-lh, polytomy resolution): kb/issues/H-timetree-mass-sizing-node-times-break-downstream-invariants.md"]
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
        .all(|node| inference.posterior[&node.key()].distribution.is_some())
    );

    let edge_lh = compute_coalescent_total_lh(
      &state.graph,
      &tc,
      &inference.coalescent_node_times()?,
      &BTreeMap::new(),
      &NoopProgress,
    )?;
    let model = CoalescentModel::new(
      &compute_lineage_counts(&state.graph, &inference.coalescent_node_times()?)?,
      &tc,
    )?;
    let node_lh = -state
      .graph
      .get_nodes()
      .map(|node| {
        let time = inference.posterior[&node.key()]
          .distribution
          .as_ref()
          .and_then(|distribution| distribution.likely_time().unwrap())
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

    pretty_assert_abs_diff_eq!(node_lh, edge_lh.value(), epsilon = 1e-10);

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
      .get_mut(&root_key)
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
      .get_mut(&root_key)
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
  #[ignore = "mass-sized node times break downstream invariants (positional log-lh, polytomy resolution): kb/issues/H-timetree-mass-sizing-node-times-break-downstream-invariants.md"]
  fn test_round_unchanged_topology_recomputes_missing_time() -> Result<(), Report> {
    let (state, context) = helpers::create_polytomy_state()?;
    let tc = Distribution::constant(10.0);
    let (mut state, _) = helpers::refine(&context, state, Some(&tc))?;
    let root_key = state.graph.get_exactly_one_root()?.key();
    state
      .time_inference
      .posterior
      .get_mut(&root_key)
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

  mod helpers {
    use super::*;

    pub(super) fn newick(state: &RoundState) -> Result<String, Report> {
      nwk_write_str(
        &state.graph,
        &state.names,
        &state.branch_lengths,
        &NwkWriteOptions::default(),
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
      let nwk_parsed = nwk_read_str("(A:0.01,B:0.01,C:0.01)root;")?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;
      let alphabet = Alphabet::new(AlphabetName::Nuc)?;
      let aln: Vec<AlignmentRecord> = read_many_fasta_str(
        indoc! {r#"
          >A
          ACGTACGTACGT
          >B
          ACGTACGTACGT
          >C
          ACGTACGTACGT
        "#},
        &alphabet,
      )?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
      let (partition, _) = MarginalReconstruction::Dense(DenseReconstruction {
        partition: PartitionMarginalDense::new(0, alphabet, &graph, &node_seq_inputs(&graph, &names, aln))?,
        gtr: jc69(JC69Params::default())?,
        node_states: BTreeMap::new(),
        edges: MarginalEdges::default(),
      })
      .marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
      let branch_model = BranchModel::Marginal(partition);

      let dates: DatesMap = btreemap! {
        "A".to_owned() => Some(DateConstraint::exact(2010.0)),
        "B".to_owned() => Some(DateConstraint::exact(2015.0)),
        "C".to_owned() => Some(DateConstraint::exact(2020.0)),
      };
      let constraints = load_date_constraints(&dates, &graph, &names, &NoopProgress)?;

      let times = likely_times(&graph, &constraints, None)?;
      let inputs = ClockInputs::from_times(&graph, &times, &BTreeMap::new());
      let (
        ClockTree {
          graph, branch_lengths, ..
        },
        clock_reroot,
      ) = estimate_clock_model_with_reroot_policy(
        ClockTree {
          graph,
          branch_lengths,
          inputs,
        },
        &BTreeSet::new(),
        &ClockVarianceParams::default(),
        Some(CLOCK_RATE),
        true,
        &BranchPointOptimizationParams::default(),
        &RerootParams::default(),
        None,
        &names,
        &NoopProgress,
      )?;
      let ClockFit {
        model: clock_model,
        points: clock_points,
      } = clock_reroot.into_clock_fit()?;
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
          no_indels: false,
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
        clock_branch_lengths: BTreeMap::new(),
        gammas,
        time_inference,
      };
      Ok((state, context(constraints)))
    }

    pub(super) fn refine(
      context: &TimetreeContext,
      state: RoundState,
      coalescent_tc: Option<&Distribution>,
    ) -> Result<(RoundState, RoundOutcome), Report> {
      let pinned_tc = Distribution::constant(ROUND_TEST_TC);
      let coalescent = CoalescentModel::new(
        &compute_lineage_counts(&state.graph, &state.time_inference.coalescent_node_times()?)?,
        coalescent_tc.unwrap_or(&pinned_tc),
      )?;
      let merger_rate =
        coalescent.branch_merger_rate_schedule(&PiecewiseConstantFn::new(array![], array![ROUND_TEST_TC]))?;
      let leaf_bad_branches = bad_leaves(&state.graph, &context.date_constraints, &BTreeSet::new());
      let outliers = BTreeSet::new();
      let params = TimetreeParams {
        clock_rate: Some(CLOCK_RATE),
        resolve_polytomies: true,
        seed: Some(ROUND_TEST_SEED),
        ..marginal_timetree_params()
      };
      let inputs = RoundInputs {
        params: &params,
        context,
        leaf_bad_branches: &leaf_bad_branches,
        outliers: &outliers,
      };

      refinement_round(
        &inputs,
        &merger_rate,
        coalescent_tc.is_some().then_some(&coalescent),
        state,
        &mut get_random_number_generator(Some(ROUND_TEST_SEED)),
        &NoopProgress,
      )
    }

    fn context(date_constraints: DateConstraints) -> TimetreeContext {
      TimetreeContext {
        time_marginal: TimeMarginalMode::Never,
        date_constraints,
        covariation_clock_params: ClockVarianceParams::default(),
        branch_params: BranchPointOptimizationParams::default(),
      }
    }
  }
}
