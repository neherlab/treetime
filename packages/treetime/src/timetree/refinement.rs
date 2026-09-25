use crate::clock::clock_model::ClockModel;
use crate::clock::clock_regression::{
  ClockRegressionPoint, ClockVarianceParams, estimate_clock_model_with_reroot_policy,
};
use crate::clock::clock_state::{ClockInputs, ClockState};
use crate::clock::date_constraints::DateConstraints;
use crate::clock::find_best_root::params::BranchPointOptimizationParams;
use crate::clock::reroot::RerootParams;
use crate::coalescent::coalescent::CoalescentModel;
use crate::partition::timetree::marginal::marginal_update_timetree;
use crate::partition::timetree::partition::PartitionTimetree;
use crate::progress::ProgressSink;
use crate::progress_info;
use crate::timetree::convergence::node_times::{NodeTimeChange, capture_node_times, measure_node_time_change};
use crate::timetree::convergence::sequence_changes::{capture_ancestral_states, count_sequence_changes};
use crate::timetree::inference::runner::{
  CLOCK_BRANCH_LENGTH_DAMPING, commit_clock_branch_lengths, run_timetree, timetree_branch_lengths,
};
use crate::timetree::inference::time_inference::{TimeInference, likely_times, unit_gammas};
use crate::timetree::optimization::polytomy::resolve::{require_internal_node_times, resolve_polytomies};
use crate::timetree::optimization::relaxed_clock::{RelaxedClockPrior, apply_relaxed_clock};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use std::collections::BTreeMap;
use treetime_graph::assign_node_names::assign_node_names;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;

pub(crate) struct Refinement<'a> {
  pub graph: &'a mut Graph,
  pub partitions: Vec<PartitionTimetree>,
  pub clock_model: &'a mut ClockModel,
  pub clock_points: &'a mut Vec<ClockRegressionPoint>,
  pub clock_params: &'a ClockVarianceParams,
  pub branch_params: &'a BranchPointOptimizationParams,
  pub merger_rate: &'a PiecewiseConstantFn,
  pub prior: Option<&'a CoalescentModel>,
  pub rng: &'a mut dyn rand::RngCore,
  pub options: &'a RefinementOptions,
  pub constraints: &'a DateConstraints,
  pub leaf_bad_branches: &'a BTreeMap<GraphNodeKey, bool>,
  pub inference: TimeInference,
  pub gammas: BTreeMap<GraphEdgeKey, f64>,
  pub clock_state: &'a mut ClockState,
  pub clock_branch_lengths: &'a mut BTreeMap<GraphEdgeKey, f64>,
  pub branch_lengths: &'a mut BTreeMap<GraphEdgeKey, Option<f64>>,
  pub names: &'a mut BTreeMap<GraphNodeKey, Option<String>>,
  pub progress: &'a dyn ProgressSink,
}

impl Refinement<'_> {
  pub(crate) fn run(mut self) -> Result<RefinementResult, Report> {
    let total_length = self.total_sequence_length();
    self.apply_relaxed_clock(total_length)?;

    let previous_times = capture_node_times(self.graph, &self.inference);
    let previous_states = capture_ancestral_states(self.graph, &self.partitions);
    let topology = self.refine_topology(total_length)?;
    self.rebuild_inference(topology.changed())?;

    commit_clock_branch_lengths(
      self.graph,
      self.clock_model.clock_rate(),
      CLOCK_BRANCH_LENGTH_DAMPING,
      self.clock_branch_lengths,
      &self.inference.node_times(),
      &self.gammas,
      self.progress,
    );

    let current_states = capture_ancestral_states(self.graph, &self.partitions);
    let time_change = measure_node_time_change(&previous_times, &capture_node_times(self.graph, &self.inference));

    self.update_clock_model()?;

    let outcome = RefinementOutcome {
      sequence_changes: count_sequence_changes(&previous_states, &current_states),
      time_change,
      topology,
    };
    Ok(RefinementResult {
      inference: self.inference,
      gammas: self.gammas,
      partitions: self.partitions,
      outcome,
    })
  }

  fn total_sequence_length(&self) -> usize {
    self
      .partitions
      .iter()
      .map(|partition| partition.get_sequence_length())
      .sum()
  }

  #[allow(
    clippy::as_conversions,
    reason = "count/index numeric cast is exact for the domain range"
  )]
  fn apply_relaxed_clock(&mut self, total_length: usize) -> Result<(), Report> {
    if self.options.relax.is_empty() {
      return Ok(());
    }
    if total_length == 0 {
      progress_info!(
        self.progress,
        "Skipping relaxed clock: no sequence data (partitions empty or zero-length)"
      );
      return Ok(());
    }

    let RelaxedClockPrior { slack, coupling } = RelaxedClockPrior::of(&self.options.relax);
    progress_info!(
      self.progress,
      "Applying relaxed clock with slack={slack}, coupling={coupling}"
    );
    self.gammas = apply_relaxed_clock(
      self.graph,
      self.branch_lengths,
      &self.options.relax,
      1.0 / total_length as f64,
      self.clock_model.clock_rate(),
      &self.inference.branches,
    )?;
    Ok(())
  }

  #[allow(
    clippy::as_conversions,
    reason = "count/index numeric cast is exact for the domain range"
  )]
  fn refine_topology(&mut self, total_length: usize) -> Result<TopologyOutcome, Report> {
    if self.options.topology == TopologyRefinement::Disabled {
      return Ok(TopologyOutcome::Unchanged);
    }

    let total_mutation_rate = self.clock_model.clock_rate() * total_length as f64;

    let mut node_times = self.inference.node_times();
    let merger_times = resolve_polytomies(
      self.graph,
      &self.partitions,
      total_mutation_rate,
      total_length,
      self.merger_rate,
      self.rng,
      self.branch_lengths,
      &node_times,
    )
    .wrap_err("Polytomy resolution failed")?;
    let resolved_nodes = merger_times.len();
    if resolved_nodes == 0 {
      return Ok(TopologyOutcome::Unchanged);
    }

    progress_info!(
      self.progress,
      "Resolved polytomies, introduced {resolved_nodes} new nodes"
    );
    *self.names = assign_node_names(std::mem::take(self.names), self.graph)?;
    node_times.extend(merger_times.into_iter().map(|(key, time)| (key, Some(time))));
    require_internal_node_times(self.graph, &node_times).wrap_err("Failed to prepare tree after topology change")?;
    self.gammas = unit_gammas(self.graph);
    let graph = &*self.graph;
    let partitions = std::mem::take(&mut self.partitions)
      .into_iter()
      .map(|partition| partition.reconcile_topology(graph))
      .collect_vec();
    self.partitions = partitions;

    commit_clock_branch_lengths(
      self.graph,
      self.clock_model.clock_rate(),
      1.0,
      self.clock_branch_lengths,
      &node_times,
      &self.gammas,
      self.progress,
    );

    Ok(TopologyOutcome::Changed { resolved_nodes })
  }

  fn rebuild_inference(&mut self, topology_changed: bool) -> Result<(), Report> {
    let run_branch_lengths = &*self.branch_lengths;
    let run_names = &*self.names;

    if !self.partitions.is_empty() {
      progress_info!(
        self.progress,
        "Updating ancestral sequences via marginal reconstruction"
      );
      let branch_lengths = timetree_branch_lengths(self.graph, run_branch_lengths, self.clock_branch_lengths);
      (self.partitions, _) =
        marginal_update_timetree(self.graph, &branch_lengths, std::mem::take(&mut self.partitions))?;
    }

    if topology_changed {
      progress_info!(
        self.progress,
        "Tree structure changed - rebuilding node-time state before coalescent inference"
      );
      self.inference = run_timetree(
        self.graph,
        self.constraints,
        self.leaf_bad_branches,
        &self.gammas,
        &self.partitions,
        run_branch_lengths,
        run_names,
        self.clock_model,
        None,
        self.options.no_indels,
        self.clock_state,
        self.progress,
      )
      .wrap_err("Coalescent-free timetree rebuild failed")?;
      if self.prior.is_none() {
        return Ok(());
      }
    } else {
      progress_info!(self.progress, "Updating node times via timetree inference");
    }

    self.inference = run_timetree(
      self.graph,
      self.constraints,
      self.leaf_bad_branches,
      &self.gammas,
      &self.partitions,
      run_branch_lengths,
      run_names,
      self.clock_model,
      self.prior,
      self.options.no_indels,
      self.clock_state,
      self.progress,
    )
    .wrap_err("Timetree inference failed")?;
    Ok(())
  }

  fn update_clock_model(&mut self) -> Result<(), Report> {
    let edge_inputs: BTreeMap<GraphEdgeKey, (Option<f64>, f64)> = self
      .inference
      .branches
      .iter()
      .map(|(key, branch)| (*key, (branch.time_length, self.gammas[key])))
      .collect();
    self.clock_state.reseed_transitional(self.graph);
    let mut clock_inputs = ClockInputs::new(self.graph);
    let times = likely_times(self.graph, self.constraints, Some(&self.inference))?;
    clock_inputs.reseed_from_times(self.graph, &times, &edge_inputs);
    let (new_clock_state, clock_reroot) = estimate_clock_model_with_reroot_policy(
      self.graph,
      &mut clock_inputs,
      std::mem::take(self.clock_state),
      self.clock_params,
      self.options.clock_rate,
      true,
      self.branch_params,
      &RerootParams::default(),
      self.branch_lengths,
      Some(self.clock_model.clock_rate()),
      self.names,
      self.progress,
    )
    .wrap_err("Failed to update clock model")?;
    *self.clock_state = new_clock_state;
    let fit = clock_reroot.into_clock_fit()?;
    *self.clock_model = fit.model;
    *self.clock_points = fit.points;
    Ok(())
  }
}

pub(crate) struct RefinementResult {
  pub inference: TimeInference,
  pub gammas: BTreeMap<GraphEdgeKey, f64>,
  pub partitions: Vec<PartitionTimetree>,
  pub outcome: RefinementOutcome,
}

pub(crate) struct RefinementOptions {
  pub relax: Vec<f64>,
  pub topology: TopologyRefinement,
  pub clock_rate: Option<f64>,
  pub no_indels: bool,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum TopologyRefinement {
  Disabled,
  Resolve,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct RefinementOutcome {
  pub sequence_changes: usize,
  pub time_change: NodeTimeChange,
  pub topology: TopologyOutcome,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
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
