use crate::clock::clock_model::ClockModel;
use crate::clock::clock_regression::{ClockParams, estimate_clock_model_with_reroot_policy};
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
use log::info;
use std::collections::BTreeMap;
use treetime_graph::assign_node_names::assign_node_names;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;

pub(crate) struct Refinement<'a> {
  pub graph: &'a mut Graph,
  /// The sequence partitions this round starts from, taken by value like `state`. A successful round
  /// returns them at their refreshed marginal reconstruction; a failed round returns an error and the
  /// partitions die with the round rather than surviving half-updated.
  pub partitions: Vec<PartitionTimetree>,
  pub clock_model: &'a mut ClockModel,
  pub clock_params: &'a ClockParams,
  pub branch_params: &'a BranchPointOptimizationParams,
  /// This round's per-branch coalescent merger-rate schedule.
  pub merger_rate: &'a PiecewiseConstantFn,
  /// The same model, as the prior imposed on node times -- `None` for a run that carries no
  /// coalescent prior, where the model exists only to give polytomy resolution a merger rate.
  pub prior: Option<&'a CoalescentModel>,
  /// Shared across refinement rounds so polytomy resolution draws from one continuous
  /// stream: re-seeding per round would correlate the sampled histories.
  pub rng: &'a mut dyn rand::RngCore,
  pub options: &'a RefinementOptions,
  /// Durable date constraints the date passes borrow. Never mutated by refinement.
  pub constraints: &'a DateConstraints,
  /// The date inference state this round starts from, taken by value. A successful round returns the
  /// refined state; a failed round returns an error without publishing any partial accepted state, so
  /// the caller keeps its pre-round value (A5).
  pub state: TimetreeState,
  /// Persistent per-node/per-edge clock state routed across the whole pipeline. The node divergence
  /// and outlier flag live here; the clock re-estimation reads them back
  /// from it, and each `run_timetree` refreshes the divergence into it.
  pub clock_state: &'a mut ClockState,
  /// Committed clock-constrained branch lengths keyed by edge, routed so the M-step damps against
  /// the previous round's value.
  pub clock_branch_lengths: &'a mut BTreeMap<GraphEdgeKey, f64>,
  /// Raw per-edge branch lengths routed across the loop. Polytomy
  /// resolution adds and merges edges and updates this map in place.
  pub branch_lengths: &'a mut BTreeMap<GraphEdgeKey, Option<f64>>,
  /// Per-node names routed across the loop. Polytomy resolution adds
  /// nodes and re-runs `assign_node_names`; this map is refreshed from that call's return so every
  /// later reader (this round's `run_timetree` and the pipeline's post-loop consumers) sees the
  /// current labels.
  pub names: &'a mut BTreeMap<GraphNodeKey, Option<String>>,
}

impl Refinement<'_> {
  pub fn run(mut self) -> Result<(TimetreeState, Vec<PartitionTimetree>, RefinementOutcome), Report> {
    let total_length = self.total_sequence_length();
    self.apply_relaxed_clock(total_length)?;

    // Node times are what the round moves, so they are the primary convergence signal. The
    // ancestral-state comparison is a hold-over from early v0, where internal node states were
    // fixed; it survives only as the fallback for a tree with no comparable dated nodes.
    let previous_times = capture_node_times(self.graph, &self.state);
    let previous_states = capture_ancestral_states(self.graph, &self.partitions);
    let topology = self.refine_topology(total_length)?;
    self.rebuild_inference(topology.changed())?;

    // Close the loop: the times just inferred become the lengths the next round's marginal
    // reconstruction propagates along. Damped, because each pass re-infers every time at once.
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

  fn refine_topology(&mut self, total_length: usize) -> Result<TopologyOutcome, Report> {
    if self.options.topology == TopologyRefinement::Disabled {
      return Ok(TopologyOutcome::Unchanged);
    }

    // Expected substitutions per unit time across the whole alignment.
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
    // Reset the value-resident edge fields for the new topology, the counterpart of the reset
    // `prepare_tree_after_topology_change` performs on the transitional fields.
    self.state.reset_date_edges_for_topology_change(self.graph);
    for partition in &mut self.partitions {
      partition.reconcile_topology(self.graph);
    }

    // Re-parenting invalidates the committed lengths, which describe a parent-child pair that no
    // longer exists. The sampled subtree dates every node it creates, so recommit from those
    // times rather than falling back to ML lengths for the reconstruction that follows. Undamped:
    // there is nothing meaningful to blend a moved edge with.
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
    // Use the routed per-edge length and per-node name maps for every pass below; the marginal
    // reconstruction and neither run_timetree call renames or re-lengths, so one read serves all.
    // Topology resolution and its `assign_node_names` ran before `rebuild_inference` and updated both
    // maps, so they reflect the current tree.
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
    // Re-read the clock inputs while preserving the value-resident divergence and outlier flag, then
    // re-estimate on the threaded state with the root kept (no reroot in the refinement loop). This
    // matches the standalone `estimate_clock_model_with_reroot` convenience, except the outlier flag
    // and the node dates come from the threaded values: the date passes have
    // refined the times on the date state since the last clock call.
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
  pub fn changed(self) -> bool {
    matches!(self, Self::Changed { .. })
  }

  pub fn resolved_nodes(self) -> usize {
    match self {
      Self::Unchanged => 0,
      Self::Changed { resolved_nodes } => resolved_nodes,
    }
  }
}
