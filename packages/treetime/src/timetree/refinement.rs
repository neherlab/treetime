use crate::clock::clock_model::ClockModel;
use crate::clock::clock_regression::{ClockVarianceParams, estimate_clock_model_with_reroot_policy};
use crate::clock::clock_state::{ClockInputs, ClockState};
use crate::clock::date_constraints::DateConstraints;
use crate::clock::find_best_root::params::BranchPointOptimizationParams;
use crate::clock::reroot::RerootParams;
use crate::coalescent::coalescent::CoalescentModel;
use crate::partition::timetree::marginal::marginal_update_timetree;
use crate::partition::timetree::partition::PartitionTimetree;
use crate::timetree::convergence::node_times::{NodeTimeChange, capture_node_times, measure_node_time_change};
use crate::timetree::convergence::sequence_changes::{capture_ancestral_states, count_sequence_changes};
use crate::timetree::inference::runner::{
  CLOCK_BRANCH_LENGTH_DAMPING, commit_clock_branch_lengths, run_timetree, timetree_branch_lengths,
};
use crate::timetree::optimization::clock_filter::propagate_bad_branches;
use crate::timetree::optimization::polytomy::{prepare_tree_after_topology_change, resolve_polytomies};
use crate::timetree::optimization::relaxed_clock::apply_relaxed_clock;
use crate::timetree::timetree_state::TimetreeState;
use eyre::{Report, WrapErr};
use itertools::Itertools;
use log::info;
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
  pub clock_params: &'a ClockVarianceParams,
  pub branch_params: &'a BranchPointOptimizationParams,
  pub merger_rate: &'a PiecewiseConstantFn,
  pub prior: Option<&'a CoalescentModel>,
  pub rng: &'a mut dyn rand::RngCore,
  pub options: &'a RefinementOptions,
  pub constraints: &'a DateConstraints,
  pub state: TimetreeState,
  pub clock_state: &'a mut ClockState,
  pub clock_branch_lengths: &'a mut BTreeMap<GraphEdgeKey, f64>,
  pub branch_lengths: &'a mut BTreeMap<GraphEdgeKey, Option<f64>>,
  pub names: &'a mut BTreeMap<GraphNodeKey, Option<String>>,
}

impl Refinement<'_> {
  pub(crate) fn run(mut self) -> Result<(TimetreeState, Vec<PartitionTimetree>, RefinementOutcome), Report> {
    let total_length = self.total_sequence_length();
    self.apply_relaxed_clock(total_length)?;

    let previous_times = capture_node_times(self.graph, &self.state);
    let previous_states = capture_ancestral_states(self.graph, &self.partitions);
    let topology = self.refine_topology(total_length)?;
    self.rebuild_inference(topology.changed())?;

    commit_clock_branch_lengths(
      self.graph,
      self.clock_model.clock_rate(),
      CLOCK_BRANCH_LENGTH_DAMPING,
      self.clock_branch_lengths,
      &self.state,
    );

    let current_states = capture_ancestral_states(self.graph, &self.partitions);
    let time_change = measure_node_time_change(&previous_times, &capture_node_times(self.graph, &self.state));

    self.update_clock_model()?;

    let outcome = RefinementOutcome {
      sequence_changes: count_sequence_changes(&previous_states, &current_states),
      time_change,
      topology,
    };
    Ok((self.state, self.partitions, outcome))
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
      info!("Skipping relaxed clock: no sequence data (partitions empty or zero-length)");
      return Ok(());
    }

    info!(
      "Applying relaxed clock with slack={}, coupling={}",
      self.options.relax.first().copied().unwrap_or(1.0),
      self.options.relax.get(1).copied().unwrap_or(1.0)
    );
    apply_relaxed_clock(
      self.graph,
      self.branch_lengths,
      &self.options.relax,
      1.0 / total_length as f64,
      self.clock_model.clock_rate(),
      &mut self.state,
    )
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

    let resolved_nodes = resolve_polytomies(
      self.graph,
      &self.partitions,
      total_mutation_rate,
      total_length,
      self.merger_rate,
      self.rng,
      self.branch_lengths,
      &mut self.state,
    )
    .wrap_err("Polytomy resolution failed")?;
    if resolved_nodes == 0 {
      return Ok(TopologyOutcome::Unchanged);
    }

    info!("Resolved polytomies, introduced {resolved_nodes} new nodes");
    *self.names = assign_node_names(std::mem::take(self.names), self.graph)?;
    propagate_bad_branches(self.graph, &mut self.state)?;
    prepare_tree_after_topology_change(self.graph, &mut self.state)
      .wrap_err("Failed to prepare tree after topology change")?;
    self.state.reset_date_edges_for_topology_change(self.graph);
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
      &self.state,
    );

    Ok(TopologyOutcome::Changed { resolved_nodes })
  }

  fn rebuild_inference(&mut self, topology_changed: bool) -> Result<(), Report> {
    let run_branch_lengths = &*self.branch_lengths;
    let run_names = &*self.names;

    if !self.partitions.is_empty() {
      info!("Updating ancestral sequences via marginal reconstruction");
      let branch_lengths = timetree_branch_lengths(self.graph, run_branch_lengths, self.clock_branch_lengths);
      (self.partitions, _) =
        marginal_update_timetree(self.graph, &branch_lengths, std::mem::take(&mut self.partitions))?;
    }

    if topology_changed {
      info!("Tree structure changed - rebuilding node-time state before coalescent inference");
      self.state = run_timetree(
        self.graph,
        self.constraints,
        &self.partitions,
        run_branch_lengths,
        run_names,
        self.clock_model,
        None,
        self.options.no_indels,
        std::mem::take(&mut self.state),
        self.clock_state,
      )
      .wrap_err("Coalescent-free timetree rebuild failed")?;
      if self.prior.is_none() {
        return Ok(());
      }
    } else {
      info!("Updating node times via timetree inference");
    }

    self.state = run_timetree(
      self.graph,
      self.constraints,
      &self.partitions,
      run_branch_lengths,
      run_names,
      self.clock_model,
      self.prior,
      self.options.no_indels,
      std::mem::take(&mut self.state),
      self.clock_state,
    )
    .wrap_err("Timetree inference failed")?;
    Ok(())
  }

  fn update_clock_model(&mut self) -> Result<(), Report> {
    let edge_inputs: BTreeMap<GraphEdgeKey, (Option<f64>, f64)> = self
      .state
      .edges
      .iter()
      .map(|(key, edge)| (*key, (edge.time_length, edge.gamma)))
      .collect();
    self.clock_state.reseed_transitional(self.graph);
    let mut clock_inputs = ClockInputs::new(self.graph);
    clock_inputs.reseed_from_times(self.graph, &self.state.likely_times(self.constraints), &edge_inputs);
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
    )
    .wrap_err("Failed to update clock model")?;
    *self.clock_state = new_clock_state;
    *self.clock_model = clock_reroot.into_clock_model()?;
    Ok(())
  }
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
