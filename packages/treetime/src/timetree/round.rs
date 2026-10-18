use crate::branch_lengths::one_mutation;
use crate::clock::clock_model::ClockModel;
use crate::clock::clock_regression::{
  ClockFit, ClockFitParams, ClockRegressionPoint, ClockTree, estimate_clock_model_with_reroot_policy,
};
use crate::clock::clock_state::ClockInputs;
use crate::clock::reroot::RerootParams;
use crate::coalescent::mean_field::CoalescentPrior;
use crate::progress::LogSink;
use crate::progress_info;
use crate::timetree::branch_model::BranchModel;
use crate::timetree::convergence::node_times::{NodeTimeChange, capture_node_times, measure_node_time_change};
use crate::timetree::convergence::sequence_changes::{capture_ancestral_states, count_sequence_changes};
use crate::timetree::inference::result::{TimeInference, likely_times};
use crate::timetree::inference::runner::{
  CLOCK_BRANCH_LENGTH_DAMPING, CLOCK_BRANCH_LENGTH_UNDAMPED, TimeInferenceInputs, blended_clock_branch_lengths,
  run_timetree, timetree_branch_lengths,
};
use crate::timetree::optimization::marginal_clock_rate::{
  MarginalRateInputs, MarginalRateStep, step_marginal_clock_rate,
};
use crate::timetree::optimization::polytomy::resolve::{
  PolytomyResolution, require_internal_node_times, resolve_polytomies,
};
use crate::timetree::optimization::relaxed_clock::{RelaxedClockPrior, apply_relaxed_clock, unit_gammas};
use crate::timetree::params::{TimetreeContext, TimetreeParams};
use eyre::{Report, WrapErr};
use rand::Rng;
use std::collections::{BTreeMap, BTreeSet};
use std::mem;
use treetime_graph::assign_node_names::assign_node_names;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
use treetime_utils::cast::ToF64Exact;

pub(crate) fn refinement_round(
  inputs: &RoundInputs<'_>,
  merger_rate: &PiecewiseConstantFn,
  prior: Option<&CoalescentPrior>,
  state: RoundState,
  rng: &mut dyn Rng,
  log: &dyn LogSink,
) -> Result<(RoundState, RoundOutcome), Report> {
  let total_length = state.branch_model.total_length();
  let state = relax_clock(&inputs.params.relax, state, total_length, log)?;

  let previous_times = capture_node_times(&state.graph, &state.time_inference);
  let previous_states = capture_ancestral_states(&state.graph, &state.branch_model);
  let (state, topology) = refine_topology(inputs, state, total_length, merger_rate, rng, log)?;
  let state = refresh_times(inputs, prior, topology.changed(), state, log)?;

  let current_states = capture_ancestral_states(&state.graph, &state.branch_model);
  let current_times = capture_node_times(&state.graph, &state.time_inference);
  let outcome = RoundOutcome {
    sequence_changes: count_sequence_changes(&previous_states, &current_states),
    time_change: measure_node_time_change(&previous_times, &current_times),
    topology,
  };
  let state = update_clock_model(inputs, prior, state, log)?;
  Ok((state, outcome))
}

pub(crate) fn infer_final_times(
  inputs: &RoundInputs<'_>,
  prior: Option<&CoalescentPrior>,
  state: &RoundState,
  log: &dyn LogSink,
) -> Result<TimeInference, Report> {
  run_timetree(&state.time_inference_inputs(inputs), prior, log).wrap_err("Final timetree inference failed")
}

pub(crate) fn final_marginal_round(
  time_inference: TimeInference,
  state: RoundState,
  log: &dyn LogSink,
) -> Result<RoundState, Report> {
  RoundState {
    time_inference,
    ..state
  }
  .blend_clock_branch_lengths(CLOCK_BRANCH_LENGTH_UNDAMPED, log)
  .marginal_update()
}

pub(crate) struct RoundOutcome {
  pub sequence_changes: usize,
  pub time_change: NodeTimeChange,
  pub topology: TopologyOutcome,
}

fn relax_clock(relax: &[f64], state: RoundState, total_length: usize, log: &dyn LogSink) -> Result<RoundState, Report> {
  if relax.is_empty() {
    return Ok(state);
  }
  if total_length == 0 {
    progress_info!(
      log,
      "Skipping relaxed clock: no sequence data (no alignment or zero sequence length)"
    );
    return Ok(state);
  }

  let RelaxedClockPrior { slack, coupling } = RelaxedClockPrior::of(relax);
  progress_info!(log, "Applying relaxed clock with slack={slack}, coupling={coupling}");
  let one_mutation = one_mutation(total_length);
  let gammas = apply_relaxed_clock(
    &state.graph,
    &state.branch_lengths,
    relax,
    one_mutation,
    state.clock_model.clock_rate(),
    &state.time_inference.branches,
  )?;
  Ok(RoundState { gammas, ..state })
}

fn refine_topology(
  inputs: &RoundInputs<'_>,
  state: RoundState,
  total_length: usize,
  merger_rate: &PiecewiseConstantFn,
  rng: &mut dyn Rng,
  log: &dyn LogSink,
) -> Result<(RoundState, TopologyOutcome), Report> {
  if !inputs.params.resolve_polytomies {
    return Ok((state, TopologyOutcome::Unchanged));
  }

  let total_mutation_rate = state.clock_model.clock_rate() * total_length.to_f64_exact();
  let mut node_times = state.time_inference.node_times();
  let PolytomyResolution {
    graph,
    branch_lengths,
    merger_times,
    removed_nodes,
  } = resolve_polytomies(
    state.graph,
    state.branch_lengths,
    &state.branch_model,
    total_mutation_rate,
    total_length,
    merger_rate,
    rng,
    &node_times,
  )
  .wrap_err("Polytomy resolution failed")?;
  if merger_times.is_empty() && removed_nodes == 0 {
    let state = RoundState {
      graph,
      branch_lengths,
      ..state
    };
    return Ok((state, TopologyOutcome::Unchanged));
  }

  let resolved_nodes = merger_times.len();
  if resolved_nodes > 0 {
    progress_info!(log, "Resolved polytomies, introduced {resolved_nodes} new nodes");
  }
  let names = assign_node_names(state.names, &graph)?.names;
  node_times.extend(merger_times.into_iter().map(|(key, time)| (key, Some(time))));
  require_internal_node_times(&graph, &node_times)
    .wrap_err("Polytomy resolution left an internal node without an inferred time")?;
  let gammas = unit_gammas(&graph);
  let branch_model = state.branch_model.reconcile_topology(&graph);
  let clock_branch_lengths = blended_clock_branch_lengths(
    &graph,
    state.clock_model.clock_rate(),
    CLOCK_BRANCH_LENGTH_UNDAMPED,
    &state.clock_branch_lengths,
    &node_times,
    &gammas,
    log,
  );

  let state = RoundState {
    graph,
    names,
    branch_model,
    branch_lengths,
    clock_branch_lengths,
    gammas,
    ..state
  };
  Ok((state, TopologyOutcome::Changed { resolved_nodes }))
}

#[derive(Clone, Copy)]
pub(crate) enum TopologyOutcome {
  Unchanged,
  Changed { resolved_nodes: usize },
}

impl TopologyOutcome {
  fn changed(self) -> bool {
    matches!(self, Self::Changed { .. })
  }

  pub(crate) fn resolved_nodes(self) -> usize {
    match self {
      Self::Unchanged => 0,
      Self::Changed { resolved_nodes } => resolved_nodes,
    }
  }
}

fn refresh_times(
  inputs: &RoundInputs<'_>,
  prior: Option<&CoalescentPrior>,
  topology_changed: bool,
  state: RoundState,
  log: &dyn LogSink,
) -> Result<RoundState, Report> {
  if state.branch_model.is_marginal() {
    progress_info!(log, "Updating ancestral sequences via marginal reconstruction");
  }
  let mut state = state.marginal_update()?;
  drop(mem::take(&mut state.time_inference));
  let time_inference = infer_times(inputs, prior, topology_changed, &state, log)?;
  Ok(
    RoundState {
      time_inference,
      ..state
    }
    .blend_clock_branch_lengths(CLOCK_BRANCH_LENGTH_DAMPING, log),
  )
}

fn infer_times(
  inputs: &RoundInputs<'_>,
  prior: Option<&CoalescentPrior>,
  topology_changed: bool,
  state: &RoundState,
  log: &dyn LogSink,
) -> Result<TimeInference, Report> {
  if topology_changed {
    progress_info!(
      log,
      "Tree structure changed - updating node times on the new topology via timetree inference"
    );
  } else {
    progress_info!(log, "Updating node times via timetree inference");
  }
  run_timetree(&state.time_inference_inputs(inputs), prior, log).wrap_err("Timetree inference failed")
}

/// Re-estimates the clock model on the round's times, keeping the root. With a coalescent prior and no user-fixed rate,
/// the rate takes one damped step toward the maximum of the marginal likelihood of the node times under the round's
/// prior, and the clock model is fitted with that rate fixed; otherwise the rate comes from the root-to-tip regression.
fn update_clock_model(
  inputs: &RoundInputs<'_>,
  prior: Option<&CoalescentPrior>,
  state: RoundState,
  log: &dyn LogSink,
) -> Result<RoundState, Report> {
  let clock_rate = match (inputs.params.clock_rate, prior) {
    (None, Some(prior)) => Some(marginal_clock_rate(inputs, prior, &state, log)?),
    (fixed, _) => fixed,
  };
  let edge_inputs: BTreeMap<GraphEdgeKey, (Option<f64>, f64)> = state
    .time_inference
    .branches
    .iter()
    .map(|(key, branch)| (*key, (branch.time_length, state.gammas[key])))
    .collect();
  let times = likely_times(&state.graph, &inputs.context.date_constraints, &state.time_inference)?;
  let clock_inputs = ClockInputs::from_times(&state.graph, &times, &edge_inputs);
  let previous_clock_rate = state.clock_model.clock_rate();
  let keep_root = true;
  let (
    ClockTree {
      graph, branch_lengths, ..
    },
    clock_reroot,
  ) = estimate_clock_model_with_reroot_policy(
    ClockTree {
      graph: state.graph,
      branch_lengths: state.branch_lengths,
      inputs: clock_inputs,
    },
    &ClockFitParams {
      outliers: inputs.outliers,
      clock_rate,
      keep_root,
      reroot_params: &RerootParams::default(),
      prev_clock_rate: Some(previous_clock_rate),
    },
    &state.names,
    log,
  )
  .wrap_err("Failed to update clock model")?;
  let ClockFit { model, points } = clock_reroot.into_clock_fit()?;
  Ok(RoundState {
    graph,
    branch_lengths,
    clock_model: model,
    clock_points: points,
    ..state
  })
}

fn marginal_clock_rate(
  inputs: &RoundInputs<'_>,
  prior: &CoalescentPrior,
  state: &RoundState,
  log: &dyn LogSink,
) -> Result<f64, Report> {
  let reference_rate = state.clock_model.clock_rate();
  let MarginalRateStep {
    clock_rate,
    target_rate,
    log_lh,
  } = step_marginal_clock_rate(&MarginalRateInputs {
    graph: &state.graph,
    date_constraints: &inputs.context.date_constraints,
    coalescent: Some(prior),
    bad_branches: &state.time_inference.bad_branches,
    branches: &state.time_inference.branches,
    reference_rate,
    max_points: inputs.params.max_grid_points,
  })
  .wrap_err("Failed to estimate the clock rate from the marginal likelihood of the node times")?;
  progress_info!(
    log,
    "Clock rate stepped toward the maximum of the marginal likelihood of the node times: {clock_rate:.6e} (was {reference_rate:.6e}, target {target_rate:.6e}), log_lh_time={log_lh:.2}"
  );
  Ok(clock_rate)
}

pub(crate) struct RoundState {
  pub graph: Graph,
  pub names: BTreeMap<GraphNodeKey, Option<String>>,
  pub branch_model: BranchModel,
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
  pub clock_model: ClockModel,
  pub clock_points: Vec<ClockRegressionPoint>,
  pub clock_branch_lengths: BTreeMap<GraphEdgeKey, f64>,
  pub gammas: BTreeMap<GraphEdgeKey, f64>,
  pub time_inference: TimeInference,
}

impl RoundState {
  pub(crate) fn time_inference_inputs<'a>(&'a self, inputs: &'a RoundInputs<'_>) -> TimeInferenceInputs<'a> {
    TimeInferenceInputs {
      graph: &self.graph,
      date_constraints: &inputs.context.date_constraints,
      leaf_bad_branches: inputs.leaf_bad_branches,
      gammas: &self.gammas,
      branch_model: &self.branch_model,
      branch_lengths: &self.branch_lengths,
      names: &self.names,
      clock_model: &self.clock_model,
      clock_rate_fixed: inputs.params.clock_rate.is_some(),
      indel_policy: inputs.params.indel_policy,
      max_grid_points: inputs.params.max_grid_points,
    }
  }

  pub(crate) fn blend_clock_branch_lengths(self, damping: f64, log: &dyn LogSink) -> Self {
    let clock_branch_lengths = blended_clock_branch_lengths(
      &self.graph,
      self.clock_model.clock_rate(),
      damping,
      &self.clock_branch_lengths,
      &self.time_inference.node_times(),
      &self.gammas,
      log,
    );
    Self {
      clock_branch_lengths,
      ..self
    }
  }

  fn marginal_update(self) -> Result<Self, Report> {
    let timetree_lengths = timetree_branch_lengths(&self.graph, &self.branch_lengths, &self.clock_branch_lengths);
    let branch_model = self.branch_model.marginal_update(&self.graph, &timetree_lengths)?;
    Ok(Self { branch_model, ..self })
  }
}

pub(crate) struct RoundInputs<'a> {
  pub params: &'a TimetreeParams,
  pub context: &'a TimetreeContext,
  pub leaf_bad_branches: &'a BTreeMap<GraphNodeKey, bool>,
  pub outliers: &'a BTreeSet<GraphNodeKey>,
}
