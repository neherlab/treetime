use crate::ancestral::marginal::update_marginal;
use crate::clock::clock_model::ClockModel;
use crate::clock::clock_regression::{ClockParams, estimate_clock_model_with_reroot};
use crate::clock::find_best_root::params::BranchPointOptimizationParams;
use crate::coalescent::coalescent::CoalescentModel;
use crate::partition::timetree::partition::{GraphTimetree, PartitionTimetreeRef};
use crate::partition::traits::{PartitionMarginalPasses, PartitionTimetreeOps};
use crate::timetree::convergence::node_times::{NodeTimeChange, capture_node_times, measure_node_time_change};
use crate::timetree::convergence::sequence_changes::{capture_ancestral_states, count_sequence_changes};
use crate::timetree::inference::runner::{CLOCK_BRANCH_LENGTH_DAMPING, commit_clock_branch_lengths, run_timetree};
use crate::timetree::optimization::clock_filter::propagate_bad_branches;
use crate::timetree::optimization::polytomy::{prepare_tree_after_topology_change, resolve_polytomies};
use crate::timetree::optimization::relaxed_clock::apply_relaxed_clock;
use eyre::{Report, WrapErr};
use log::info;
use treetime_graph::assign_node_names::assign_node_names;
use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;

pub(crate) struct Refinement<'a> {
  pub graph: &'a mut GraphTimetree,
  pub partitions: &'a [PartitionTimetreeRef],
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
}

impl Refinement<'_> {
  pub fn run(mut self) -> Result<RefinementOutcome, Report> {
    let total_length = self.total_sequence_length();
    self.apply_relaxed_clock(total_length)?;

    // Node times are what the round moves, so they are the primary convergence signal. The
    // ancestral-state comparison is a hold-over from early v0, where internal node states were
    // fixed; it survives only as the fallback for a tree with no comparable dated nodes.
    let previous_times = capture_node_times(self.graph);
    let previous_states = capture_ancestral_states(self.graph, self.partitions);
    let topology = self.refine_topology(total_length)?;
    self.rebuild_inference(topology.changed())?;

    // Close the loop: the times just inferred become the lengths the next round's marginal
    // reconstruction propagates along. Damped, because each pass re-infers every time at once.
    commit_clock_branch_lengths(self.graph, self.clock_model.clock_rate(), CLOCK_BRANCH_LENGTH_DAMPING);

    let current_states = capture_ancestral_states(self.graph, self.partitions);
    let time_change = measure_node_time_change(&previous_times, &capture_node_times(self.graph));

    self.update_clock_model()?;

    Ok(RefinementOutcome {
      sequence_changes: count_sequence_changes(&previous_states, &current_states),
      time_change,
      topology,
    })
  }

  fn total_sequence_length(&self) -> usize {
    self
      .partitions
      .iter()
      .map(|partition| partition.read_arc().get_sequence_length())
      .sum()
  }

  fn apply_relaxed_clock(&self, total_length: usize) -> Result<(), Report> {
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
      &self.options.relax,
      1.0 / total_length as f64,
      self.clock_model.clock_rate(),
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
      self.partitions,
      total_mutation_rate,
      total_length,
      self.merger_rate,
      self.rng,
    )
    .wrap_err("Polytomy resolution failed")?;
    if resolved_nodes == 0 {
      return Ok(TopologyOutcome::Unchanged);
    }

    info!("Resolved polytomies, introduced {resolved_nodes} new nodes");
    assign_node_names(self.graph)?;
    propagate_bad_branches(self.graph)?;
    prepare_tree_after_topology_change(self.graph).wrap_err("Failed to prepare tree after topology change")?;
    for partition in self.partitions {
      partition.write_arc().reconcile_topology(self.graph);
    }

    // Re-parenting invalidates the committed lengths, which describe a parent-child pair that no
    // longer exists. The sampled subtree dates every node it creates, so recommit from those
    // times rather than falling back to ML lengths for the reconstruction that follows. Undamped:
    // there is nothing meaningful to blend a moved edge with.
    commit_clock_branch_lengths(self.graph, self.clock_model.clock_rate(), 1.0);

    Ok(TopologyOutcome::Changed { resolved_nodes })
  }

  fn rebuild_inference(&mut self, topology_changed: bool) -> Result<(), Report> {
    if !self.partitions.is_empty() {
      info!("Updating ancestral sequences via marginal reconstruction");
      update_marginal(self.graph, self.partitions)?;
    }

    if topology_changed {
      info!("Tree structure changed - rebuilding node-time state before coalescent inference");
      run_timetree(
        self.graph,
        self.partitions,
        self.clock_model,
        None,
        self.options.no_indels,
      )
      .wrap_err("Coalescent-free timetree rebuild failed")?;
      if self.prior.is_none() {
        return Ok(());
      }
    } else {
      info!("Updating node times via timetree inference");
    }

    run_timetree(
      self.graph,
      self.partitions,
      self.clock_model,
      self.prior,
      self.options.no_indels,
    )
    .wrap_err("Timetree inference failed")
  }

  fn update_clock_model(&mut self) -> Result<(), Report> {
    *self.clock_model = estimate_clock_model_with_reroot(
      self.graph,
      self.clock_params,
      self.options.clock_rate,
      true,
      self.branch_params,
      Some(self.clock_model.clock_rate()),
    )
    .wrap_err("Failed to update clock model")?;
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
