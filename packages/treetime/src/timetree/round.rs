use crate::clock::clock_model::ClockModel;
use crate::clock::clock_regression::{
  ClockFit, ClockRegressionPoint, ClockTree, ClockVarianceParams, estimate_clock_model_with_reroot_policy,
};
use crate::clock::clock_state::ClockInputs;
use crate::clock::date_constraints::DateConstraints;
use crate::clock::find_best_root::params::BranchPointOptimizationParams;
use crate::clock::reroot::RerootParams;
use crate::coalescent::coalescent::CoalescentModel;
use crate::progress::ProgressSink;
use crate::progress_info;
use crate::timetree::branch_model::BranchModel;
use crate::timetree::convergence::node_times::{NodeTimeChange, capture_node_times, measure_node_time_change};
use crate::timetree::convergence::sequence_changes::{capture_ancestral_states, count_sequence_changes};
use crate::timetree::inference::runner::{
  CLOCK_BRANCH_LENGTH_DAMPING, CLOCK_BRANCH_LENGTH_UNDAMPED, blended_clock_branch_lengths, run_timetree,
  timetree_branch_lengths,
};
use crate::timetree::inference::time_inference::{TimeInference, likely_times, unit_gammas};
use crate::timetree::optimization::polytomy::resolve::{
  PolytomyResolution, require_internal_node_times, resolve_polytomies,
};
use crate::timetree::optimization::relaxed_clock::{RelaxedClockPrior, apply_relaxed_clock};
use crate::timetree::pipeline::{TimetreeContext, TimetreeParams};
use eyre::{Report, WrapErr};
use rand::RngCore;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::assign_node_names::assign_node_names;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;

pub(crate) fn refinement_round(
  inputs: &RoundInputs<'_>,
  merger_rate: &PiecewiseConstantFn,
  prior: Option<&CoalescentModel>,
  state: RoundState,
  rng: &mut dyn RngCore,
  progress: &dyn ProgressSink,
) -> Result<(RoundState, RoundOutcome), Report> {
  let total_length = state.branch_model.sequence_length();
  let state = relax_clock(inputs.relax, state, total_length, progress)?;

  let previous_times = capture_node_times(&state.graph, &state.time_inference);
  let previous_states = capture_ancestral_states(&state.graph, &state.branch_model);
  let (state, topology) = refine_topology(inputs.topology, state, total_length, merger_rate, rng, progress)?;
  let state = refresh_times(inputs, prior, topology.changed(), state, progress)?;

  let current_states = capture_ancestral_states(&state.graph, &state.branch_model);
  let current_times = capture_node_times(&state.graph, &state.time_inference);
  let outcome = RoundOutcome {
    sequence_changes: count_sequence_changes(&previous_states, &current_states),
    time_change: measure_node_time_change(&previous_times, &current_times),
    topology,
  };
  let state = update_clock_model(inputs, state, progress)?;
  Ok((state, outcome))
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

pub(crate) struct RoundInputs<'a> {
  constraints: &'a DateConstraints,
  leaf_bad_branches: &'a BTreeMap<GraphNodeKey, bool>,
  outliers: &'a BTreeSet<GraphNodeKey>,
  clock_params: &'a ClockVarianceParams,
  branch_params: &'a BranchPointOptimizationParams,
  relax: &'a [f64],
  topology: TopologyRefinement,
  clock_rate: Option<f64>,
  no_indels: bool,
}

impl<'a> RoundInputs<'a> {
  pub(crate) fn new(
    params: &'a TimetreeParams,
    context: &'a TimetreeContext,
    leaf_bad_branches: &'a BTreeMap<GraphNodeKey, bool>,
    outliers: &'a BTreeSet<GraphNodeKey>,
  ) -> Self {
    Self {
      constraints: &context.date_constraints,
      leaf_bad_branches,
      outliers,
      clock_params: &context.covariation_clock_params,
      branch_params: &context.branch_params,
      relax: &params.relax,
      topology: if params.resolve_polytomies {
        TopologyRefinement::Resolve
      } else {
        TopologyRefinement::Disabled
      },
      clock_rate: params.clock_rate,
      no_indels: params.no_indels,
    }
  }
}

#[derive(Clone, Copy)]
enum TopologyRefinement {
  Disabled,
  Resolve,
}

pub(crate) struct RoundOutcome {
  pub sequence_changes: usize,
  pub time_change: NodeTimeChange,
  pub topology: TopologyOutcome,
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

#[allow(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
fn relax_clock(
  relax: &[f64],
  state: RoundState,
  total_length: usize,
  progress: &dyn ProgressSink,
) -> Result<RoundState, Report> {
  if relax.is_empty() {
    return Ok(state);
  }
  if total_length == 0 {
    progress_info!(
      progress,
      "Skipping relaxed clock: no sequence data (partitions empty or zero-length)"
    );
    return Ok(state);
  }

  let RelaxedClockPrior { slack, coupling } = RelaxedClockPrior::of(relax);
  progress_info!(
    progress,
    "Applying relaxed clock with slack={slack}, coupling={coupling}"
  );
  let gammas = apply_relaxed_clock(
    &state.graph,
    &state.branch_lengths,
    relax,
    1.0 / total_length as f64,
    state.clock_model.clock_rate(),
    &state.time_inference.branches,
  )?;
  Ok(RoundState { gammas, ..state })
}

#[allow(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
fn refine_topology(
  topology: TopologyRefinement,
  state: RoundState,
  total_length: usize,
  merger_rate: &PiecewiseConstantFn,
  rng: &mut dyn RngCore,
  progress: &dyn ProgressSink,
) -> Result<(RoundState, TopologyOutcome), Report> {
  if matches!(topology, TopologyRefinement::Disabled) {
    return Ok((state, TopologyOutcome::Unchanged));
  }

  let total_mutation_rate = state.clock_model.clock_rate() * total_length as f64;
  let mut node_times = state.time_inference.node_times();
  let PolytomyResolution {
    graph,
    branch_lengths,
    merger_times,
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
  let resolved_nodes = merger_times.len();
  if resolved_nodes == 0 {
    let state = RoundState {
      graph,
      branch_lengths,
      ..state
    };
    return Ok((state, TopologyOutcome::Unchanged));
  }

  progress_info!(progress, "Resolved polytomies, introduced {resolved_nodes} new nodes");
  let names = assign_node_names(state.names, &graph)?;
  node_times.extend(merger_times.into_iter().map(|(key, time)| (key, Some(time))));
  require_internal_node_times(&graph, &node_times).wrap_err("Failed to prepare tree after topology change")?;
  let gammas = unit_gammas(&graph);
  let branch_model = state.branch_model.reconcile_topology(&graph);
  let clock_branch_lengths = blended_clock_branch_lengths(
    &graph,
    state.clock_model.clock_rate(),
    CLOCK_BRANCH_LENGTH_UNDAMPED,
    &state.clock_branch_lengths,
    &node_times,
    &gammas,
    progress,
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

fn refresh_times(
  inputs: &RoundInputs<'_>,
  prior: Option<&CoalescentModel>,
  topology_changed: bool,
  state: RoundState,
  progress: &dyn ProgressSink,
) -> Result<RoundState, Report> {
  if state.branch_model.is_marginal() {
    progress_info!(progress, "Updating ancestral sequences via marginal reconstruction");
  }
  let timetree_lengths = timetree_branch_lengths(&state.graph, &state.branch_lengths, &state.clock_branch_lengths);
  let branch_model = state.branch_model.marginal_update(&state.graph, &timetree_lengths)?;
  let state = RoundState { branch_model, ..state };
  let time_inference = infer_times(inputs, prior, topology_changed, &state, progress)?;

  let clock_branch_lengths = blended_clock_branch_lengths(
    &state.graph,
    state.clock_model.clock_rate(),
    CLOCK_BRANCH_LENGTH_DAMPING,
    &state.clock_branch_lengths,
    &time_inference.node_times(),
    &state.gammas,
    progress,
  );
  Ok(RoundState {
    clock_branch_lengths,
    time_inference,
    ..state
  })
}

fn infer_times(
  inputs: &RoundInputs<'_>,
  prior: Option<&CoalescentModel>,
  topology_changed: bool,
  state: &RoundState,
  progress: &dyn ProgressSink,
) -> Result<TimeInference, Report> {
  let run = |prior: Option<&CoalescentModel>| {
    run_timetree(
      &state.graph,
      inputs.constraints,
      inputs.leaf_bad_branches,
      &state.gammas,
      &state.branch_model,
      &state.branch_lengths,
      &state.names,
      &state.clock_model,
      prior,
      inputs.no_indels,
      progress,
    )
  };

  if topology_changed {
    progress_info!(
      progress,
      "Tree structure changed - rebuilding node-time state before coalescent inference"
    );
    let time_inference = run(None).wrap_err("Coalescent-free timetree rebuild failed")?;
    if prior.is_none() {
      return Ok(time_inference);
    }
  } else {
    progress_info!(progress, "Updating node times via timetree inference");
  }

  run(prior).wrap_err("Timetree inference failed")
}

fn update_clock_model(
  inputs: &RoundInputs<'_>,
  state: RoundState,
  progress: &dyn ProgressSink,
) -> Result<RoundState, Report> {
  let edge_inputs: BTreeMap<GraphEdgeKey, (Option<f64>, f64)> = state
    .time_inference
    .branches
    .iter()
    .map(|(key, branch)| (*key, (branch.time_length, state.gammas[key])))
    .collect();
  let times = likely_times(&state.graph, inputs.constraints, Some(&state.time_inference))?;
  let clock_inputs = ClockInputs::from_times(&state.graph, &times, &edge_inputs);
  let previous_clock_rate = state.clock_model.clock_rate();
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
    inputs.outliers,
    inputs.clock_params,
    inputs.clock_rate,
    true,
    inputs.branch_params,
    &RerootParams::default(),
    Some(previous_clock_rate),
    &state.names,
    progress,
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
