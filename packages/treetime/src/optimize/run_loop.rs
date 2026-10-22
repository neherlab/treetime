use crate::branch_lengths::branch_lengths_or_zero;
use crate::optimize::branch_length::invalid_branch_length_descriptions;
use crate::optimize::dispatch::{guess_initial_branch_lengths, optimize_edge_lengths};
use crate::optimize::gather::InitialGuessEvidence;
use crate::optimize::indel::EdgeIndels;
use crate::optimize::iteration::apply_damping;
use crate::optimize::params::ExistingBranchLengths;
use crate::optimize::params::{InitialGuessMode, OptimizeLoopSettings, TopologyOps};
use crate::optimize::partition_rates::update_partition_rates;
use crate::optimize::solver_outcome::ConvergenceSummary;
use crate::optimize::topology::collapse::collapse_edge;
use crate::optimize::topology::resolve_polytomy::resolve_polytomies;
use crate::partition::set::PartitionSet;
use crate::progress::LogSink;
use crate::{progress_info, progress_warn};
use eyre::Report;
use itertools::Itertools;
use log::debug;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node_names::NodeNames;
use treetime_primitives::LogLh;
use treetime_utils::fmt::float::float_to_significant_digits;
use treetime_utils::{make_error, o};

pub(crate) fn run_optimize_loop(
  graph: &mut Graph,
  mut partitions: PartitionSet,
  settings: &OptimizeLoopSettings,
  branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
  names: NodeNames,
  log: &dyn LogSink,
) -> Result<OptimizeLoopResult, Report> {
  let OptimizeLoopSettings {
    max_iter,
    dp,
    indel_policy,
    ..
  } = *settings;
  let mut state = LoopState {
    branch_lengths,
    names,
    convergence: ConvergenceSummary::default(),
  };
  let mut indels = indel_policy.edge_indels(graph, &partitions, &state.branch_lengths);

  let mut lh_history: Vec<LogLh> = Vec::with_capacity(max_iter + 1);
  let mut stopped_at: Option<(usize, ConvergenceReason)> = None;
  let mut best: Option<BestIteration> = None;

  for i in 0..max_iter {
    let iteration = compute_iteration(graph, &state.branch_lengths, partitions, &indels)?;
    partitions = iteration.partitions;
    let lh = iteration.likelihood.total_lh;
    lh_history.push(lh);
    log_iteration_likelihood(i + 1, &iteration.likelihood);

    if !lh.value().is_finite() {
      stopped_at = Some((i, ConvergenceReason::NumericalFailure));
      break;
    }

    record_if_best(&mut best, i, lh, &state.branch_lengths, &partitions, &indels);
    let best_lh = best.as_ref().map_or(lh, |best| best.log_lh);
    if let Some(reason) = stop_reason(&lh_history, best_lh, dp) {
      stopped_at = Some((i, reason));
      break;
    }

    let cleanup = update_branch_lengths(graph, partitions, i, settings, &indels, &mut state)?;
    if cleanup.topology_changed {
      best = None;
    }
    partitions = update_partition_rates(
      graph,
      cleanup.partitions,
      settings.partition_rates,
      &mut state.branch_lengths,
      Some(&mut indels),
      log,
    )?;
  }

  if stopped_at.is_none() && max_iter > 0 {
    let iteration = compute_iteration(graph, &state.branch_lengths, partitions, &indels)?;
    partitions = iteration.partitions;
    let lh = iteration.likelihood.total_lh;
    lh_history.push(lh);
    log_iteration_likelihood(max_iter + 1, &iteration.likelihood);
    if lh.value().is_finite() {
      record_if_best(&mut best, max_iter, lh, &state.branch_lengths, &partitions, &indels);
    }
  }

  log_stop(log, stopped_at, max_iter, dp);

  if let Some(best) = best
    && best.iteration + 1 != lh_history.len()
  {
    partitions = best.restore(graph, partitions, &mut state.branch_lengths, &mut indels)?;
  }

  Ok(OptimizeLoopResult {
    indel_rates: indels.partition_rates(partitions.len()),
    partitions,
    branch_lengths: state.branch_lengths,
    names: state.names,
    lh_history,
    stopped_at,
    convergence: state.convergence,
  })
}

#[derive(Clone, Debug)]
pub(crate) struct OptimizeLoopResult {
  pub partitions: PartitionSet,

  pub indel_rates: Vec<Option<f64>>,

  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,

  pub names: NodeNames,

  #[cfg_attr(not(test), expect(dead_code, reason = "read only in test code"))]
  pub lh_history: Vec<LogLh>,

  #[cfg_attr(not(test), expect(dead_code, reason = "read only in test code"))]
  pub stopped_at: Option<(usize, ConvergenceReason)>,

  pub convergence: ConvergenceSummary,
}

pub(super) fn record_if_best(
  best: &mut Option<BestIteration>,
  iteration: usize,
  log_lh: LogLh,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  partitions: &PartitionSet,
  indels: &EdgeIndels,
) {
  if best.as_ref().is_none_or(|best| log_lh > best.log_lh) {
    *best = Some(BestIteration {
      iteration,
      log_lh,
      branch_lengths: branch_lengths.clone(),
      rates: partitions.rates(),
      indels: indels.clone(),
    });
  }
}

pub(super) struct BestIteration {
  iteration: usize,
  log_lh: LogLh,
  branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
  rates: Vec<f64>,
  indels: EdgeIndels,
}

impl BestIteration {
  pub(super) fn restore(
    self,
    graph: &Graph,
    mut partitions: PartitionSet,
    branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
    indels: &mut EdgeIndels,
  ) -> Result<PartitionSet, Report> {
    debug!(
      "Returning the branch lengths and rates of iteration {} with the highest likelihood {}",
      self.iteration + 1,
      float_to_significant_digits(self.log_lh.value(), 7)
    );
    *branch_lengths = self.branch_lengths;
    partitions.set_rates(&self.rates)?;
    *indels = self.indels;
    Ok(
      partitions
        .marginal_update(graph, &branch_lengths_or_zero(branch_lengths))?
        .0,
    )
  }
}

fn log_stop(log: &dyn LogSink, stopped_at: Option<(usize, ConvergenceReason)>, max_iter: usize, dp: f64) {
  match stopped_at {
    Some((i, reason)) => {
      let cause = match reason {
        ConvergenceReason::Converged => format!("the log-likelihood changed by less than {dp}"),
        ConvergenceReason::Oscillating => {
          format!("the log-likelihood is within {dp} of the iteration before the previous one")
        },
        ConvergenceReason::Worsened => o!("the log-likelihood decreased"),
        ConvergenceReason::NumericalFailure => o!("the log-likelihood is not finite"),
      };
      progress_info!(log, "Optimization loop stopped at iteration {}: {cause}", i + 1);
    },
    None if max_iter == 0 => progress_info!(log, "Optimization loop skipped: the iteration limit is 0"),
    None => progress_info!(log, "Optimization loop reached the iteration limit of {max_iter}"),
  }
}

fn stop_reason(lh_history: &[LogLh], best_lh: LogLh, dp: f64) -> Option<ConvergenceReason> {
  let iteration = lh_history.len().checked_sub(1)?;
  let earlier = |steps: usize| {
    iteration
      .checked_sub(steps)
      .map_or(LogLh::IMPOSSIBLE, |i| lh_history[i])
  };
  let lh = lh_history[iteration];
  let lh_prev = earlier(1);
  let lh_prev_prev = earlier(2);
  if (lh - lh_prev).abs() < dp.abs() {
    Some(ConvergenceReason::Converged)
  } else if iteration >= 2 && (lh - lh_prev_prev).abs() < dp.abs() {
    Some(ConvergenceReason::Oscillating)
  } else if iteration >= 2 && lh < lh_prev && lh_prev >= best_lh {
    Some(ConvergenceReason::Worsened)
  } else {
    None
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ConvergenceReason {
  Converged,
  Oscillating,
  Worsened,
  NumericalFailure,
}

fn update_branch_lengths(
  graph: &mut Graph,
  partitions: PartitionSet,
  iteration: usize,
  settings: &OptimizeLoopSettings,
  indels: &EdgeIndels,
  state: &mut LoopState,
) -> Result<TopologyCleanup, Report> {
  let LoopState {
    branch_lengths,
    names,
    convergence,
  } = state;
  let old_branch_lengths = branch_lengths.clone();
  let unconverged = optimize_edge_lengths(
    graph,
    partitions.total_length(),
    &partitions,
    &indels.recount(graph, &partitions),
    settings.opt_method,
    branch_lengths,
  )?;
  convergence.record(&unconverged, graph, names);

  let zero_optimal_edges = if settings.topology_ops.collapse_short_branches {
    find_zero_optimal_internal_edges(graph, &partitions, branch_lengths)
  } else {
    vec![]
  };

  apply_damping(branch_lengths, &old_branch_lengths, settings.damping, iteration);

  prune_and_merge_in_loop(
    graph,
    partitions,
    &zero_optimal_edges,
    settings.topology_ops,
    branch_lengths,
    names,
  )
}

struct LoopState {
  branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
  names: NodeNames,
  convergence: ConvergenceSummary,
}

fn log_iteration_likelihood(iteration_number: usize, likelihood: &OptimizeIterationLikelihood) {
  debug!(
    "Iteration {iteration_number}: likelihood {} (sparse: {}, dense: {}, indel: {}, rate: {})",
    float_to_significant_digits(likelihood.total_lh.value(), 7),
    float_to_significant_digits(likelihood.sparse_lh.value(), 7),
    float_to_significant_digits(likelihood.dense_lh.value(), 7),
    float_to_significant_digits(likelihood.indel_lh.value(), 7),
    float_to_significant_digits(likelihood.indel_rate, 7)
  );
}

fn compute_iteration(
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  partitions: PartitionSet,
  indels: &EdgeIndels,
) -> Result<OptimizeIteration, Report> {
  let (partitions, marginal_lh) = partitions.marginal_update(graph, &branch_lengths_or_zero(branch_lengths))?;
  let (sparse_lh, dense_lh) = match partitions {
    PartitionSet::Sparse(_) => (marginal_lh, LogLh::ZERO),
    PartitionSet::Dense(_) => (LogLh::ZERO, marginal_lh),
  };
  let indel_lh = indels.recount(graph, &partitions).total_log_lh(graph, branch_lengths)?;
  let total_lh = sparse_lh + dense_lh + indel_lh;

  Ok(OptimizeIteration {
    partitions,
    likelihood: OptimizeIterationLikelihood {
      sparse_lh,
      dense_lh,
      indel_lh,
      total_lh,
      indel_rate: indels.rate(),
    },
  })
}

struct OptimizeIteration {
  partitions: PartitionSet,
  likelihood: OptimizeIterationLikelihood,
}

#[derive(Clone, Copy, Debug, Default)]
struct OptimizeIterationLikelihood {
  sparse_lh: LogLh,
  dense_lh: LogLh,
  indel_lh: LogLh,
  total_lh: LogLh,
  indel_rate: f64,
}

pub(crate) fn find_zero_optimal_internal_edges(
  graph: &Graph,
  partitions: &PartitionSet,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Vec<GraphEdgeKey> {
  graph
    .get_edges()
    .filter_map(|edge| {
      let bl = branch_lengths[&edge.key()].unwrap_or(f64::NAN);
      if bl != 0.0 || graph.is_leaf(edge.target()) {
        return None;
      }
      let edge_key = edge.key();
      let has_mutations = partitions
        .iter()
        .any(|partition| partition.edge_has_fitch_mutations(edge_key));
      (!has_mutations).then_some(edge_key)
    })
    .collect_vec()
}

pub(crate) fn prune_and_merge_in_loop(
  graph: &mut Graph,
  mut partitions: PartitionSet,
  zero_optimal_edges: &[GraphEdgeKey],
  topology_ops: TopologyOps,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  names: &mut NodeNames,
) -> Result<TopologyCleanup, Report> {
  let sparse = partitions.sparse_mut();

  let mut topology_changed = false;

  if !zero_optimal_edges.is_empty() {
    for &edge_key in zero_optimal_edges {
      if branch_lengths.contains_key(&edge_key) {
        branch_lengths.insert(edge_key, Some(0.0));
      }
    }

    let mut collapsed = 0_usize;
    for &edge_key in zero_optimal_edges {
      if graph.get_edge(edge_key).is_none() {
        continue;
      }
      collapse_edge(graph, sparse, edge_key, branch_lengths)?;
      collapsed += 1;
    }

    if collapsed > 0 {
      debug!("Collapsed {collapsed} zero-optimal internal edges");
      topology_changed = true;
    }
  }

  if resolve_polytomies(graph, names, sparse, topology_ops, branch_lengths)? > 0 {
    topology_changed = true;
  }

  if !topology_changed {
    return Ok(TopologyCleanup {
      partitions,
      topology_changed,
    });
  }

  graph.build()?;

  Ok(TopologyCleanup {
    partitions: partitions.reconcile_node_states(graph),
    topology_changed,
  })
}

pub(crate) struct TopologyCleanup {
  pub partitions: PartitionSet,
  pub topology_changed: bool,
}

pub(crate) fn apply_initial_guess_mode(
  graph: &Graph,
  indels: &EdgeIndels,
  evidence: &InitialGuessEvidence,
  mode: InitialGuessMode,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  names: &NodeNames,
  log: &dyn LogSink,
) -> Result<(), Report> {
  let invalid_branch_lengths = invalid_branch_length_descriptions(graph, branch_lengths, names)?;
  if let Some(message) = invalid_branch_length_warning(&invalid_branch_lengths) {
    progress_warn!(log, "{message}");
  }

  match mode {
    InitialGuessMode::Auto => {
      guess_initial_branch_lengths(graph, indels, evidence, ExistingBranchLengths::Keep, branch_lengths)
    },
    InitialGuessMode::Always => guess_initial_branch_lengths(
      graph,
      indels,
      evidence,
      ExistingBranchLengths::Overwrite,
      branch_lengths,
    ),
    InitialGuessMode::Never => {
      if !invalid_branch_lengths.is_empty() {
        return make_error!(
          "--initial-branch-lengths=never requires every edge to have a finite, non-negative branch length. \
           Invalid edges:\n  {}\n\
           Use 'auto' to fill in missing values or 'always' to recompute all branch lengths",
          invalid_branch_lengths.join("\n  ")
        );
      }
      if indels.has_indel_edge_with_zero_length(graph, branch_lengths) {
        return make_error!(
          "--initial-branch-lengths=never requires non-zero branch lengths on edges that carry indels, \
           but some indel-bearing edges have branch length zero. \
           The Poisson indel log-likelihood diverges at t=0 when k>0, so zero is not a valid input there. \
           Use 'auto' or 'always' to bootstrap these edges, or provide positive branch lengths"
        );
      }
      Ok(())
    },
  }
}

pub(super) fn invalid_branch_length_warning(invalid_branch_lengths: &[String]) -> Option<String> {
  (!invalid_branch_lengths.is_empty()).then(|| {
    format!(
      "Input tree contains invalid branch lengths (expected finite values >= 0):\n  {}\n\
       Use --initial-branch-lengths=auto to replace invalid values or \
       --initial-branch-lengths=always to recompute all branch lengths",
      invalid_branch_lengths.join("\n  ")
    )
  })
}
