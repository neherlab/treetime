use crate::branch_lengths::branch_lengths_or_zero;
use crate::branch_lengths::invalid_branch_length_descriptions;
use crate::optimize::dispatch::{guess_initial_branch_lengths, optimize_edge_lengths};
use crate::optimize::gather::InitialGuessEvidence;
use crate::optimize::indel::EdgeIndels;
use crate::optimize::iteration::apply_damping;
use crate::optimize::loop_history::{ConvergenceReason, TopologyHistory};
use crate::optimize::params::ExistingBranchLengths;
use crate::optimize::params::{InitialGuessMode, OptimizeLoopSettings, TopologyOps};
use crate::optimize::partition_rates::update_partition_rates;
use crate::optimize::solver_outcome::ConvergenceSummary;
use crate::optimize::topology::collapse::{collapse_edge, collapse_short_internal_edges};
use crate::optimize::topology::merge_shared_mutations::MergeLedger;
use crate::optimize::topology::resolve_polytomy::resolve_polytomies;
use crate::partition::marginal::sparse::alignment_indels::sparse_alignment_groups;
use crate::partition::set::PartitionSet;
use crate::progress::LogSink;
use crate::{progress_info, progress_warn};
use eyre::Report;
use itertools::Itertools;
use log::debug;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::EdgeMap;
use treetime_graph::node_names::NodeNames;
use treetime_primitives::LogLh;
use treetime_utils::fmt::float::float_to_significant_digits;
use treetime_utils::{make_error, o};

pub(crate) fn run_optimize_loop(
  graph: &mut Graph,
  mut partitions: PartitionSet,
  settings: &OptimizeLoopSettings,
  branch_lengths: EdgeMap<Option<f64>>,
  names: NodeNames,
  log: &dyn LogSink,
) -> Result<OptimizeLoopResult, Report> {
  let OptimizeLoopSettings {
    max_iter,
    dp,
    indel_policy,
    topology_ops,
    ..
  } = *settings;
  let mut state = LoopState {
    indels: indel_policy.edge_indels(graph, &partitions, &branch_lengths)?,
    branch_lengths,
    names,
    convergence: ConvergenceSummary::default(),
    ledger: MergeLedger::new(indel_policy),
    lh_history: Vec::with_capacity(max_iter + 2),
    history: TopologyHistory::new(dp),
  };

  let mut stopped_at = None;
  let mut topology_changed = false;
  for i in 0..max_iter {
    let stop;
    (partitions, stop) = evaluate(graph, partitions, &mut state)?;
    if let Some(reason) = stop {
      stopped_at = Some((i, reason));
      break;
    }
    (partitions, topology_changed) = optimize_iteration(graph, partitions, i, settings, topology_ops, &mut state, log)?;
  }

  if stopped_at.is_none() && max_iter > 0 {
    if topology_changed {
      (partitions, _) = evaluate(graph, partitions, &mut state)?;
      (partitions, _) = optimize_iteration(
        graph,
        partitions,
        max_iter,
        settings,
        TopologyOps::NONE,
        &mut state,
        log,
      )?;
    }
    (partitions, _) = evaluate(graph, partitions, &mut state)?;
  }
  log_stop(log, stopped_at, max_iter, dp);

  let LoopState {
    mut branch_lengths,
    names,
    convergence,
    mut indels,
    lh_history,
    history,
    ..
  } = state;
  if let Some(best) = history.into_best()
    && best.iteration() + 1 != lh_history.len()
  {
    partitions = best.restore(graph, partitions, &mut branch_lengths, &mut indels)?;
  }

  if topology_ops.collapse_short_branches && max_iter > 0 {
    let current_indels = indels.recount(graph, &partitions)?;
    partitions = collapse_short_internal_edges(graph, partitions, &current_indels, &mut branch_lengths)?;
  }

  Ok(OptimizeLoopResult {
    indel_rates: indels.alignment_rates(),
    partitions,
    branch_lengths,
    names,
    lh_history,
    stopped_at,
    convergence,
  })
}

#[derive(Clone, Debug)]
pub(crate) struct OptimizeLoopResult {
  pub partitions: PartitionSet,

  pub indel_rates: BTreeMap<usize, f64>,

  pub branch_lengths: EdgeMap<Option<f64>>,

  pub names: NodeNames,

  #[cfg_attr(not(test), expect(dead_code, reason = "read only in test code"))]
  pub lh_history: Vec<LogLh>,

  #[cfg_attr(not(test), expect(dead_code, reason = "read only in test code"))]
  pub stopped_at: Option<(usize, ConvergenceReason)>,

  pub convergence: ConvergenceSummary,
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

fn evaluate(
  graph: &Graph,
  partitions: PartitionSet,
  state: &mut LoopState,
) -> Result<(PartitionSet, Option<ConvergenceReason>), Report> {
  let iteration = compute_iteration(graph, &state.branch_lengths, partitions, &state.indels)?;
  let lh = iteration.likelihood.total_lh;
  let index = state.lh_history.len();
  state.lh_history.push(lh);
  log_iteration_likelihood(index + 1, &iteration.likelihood);
  let stop = state
    .history
    .record(index, lh, &state.branch_lengths, &iteration.partitions, &state.indels);
  Ok((iteration.partitions, stop))
}

fn optimize_iteration(
  graph: &mut Graph,
  partitions: PartitionSet,
  iteration: usize,
  settings: &OptimizeLoopSettings,
  topology_ops: TopologyOps,
  state: &mut LoopState,
  log: &dyn LogSink,
) -> Result<(PartitionSet, bool), Report> {
  let old_branch_lengths = state.branch_lengths.clone();
  let indels = state.indels.recount(graph, &partitions)?;
  let unconverged = optimize_edge_lengths(
    graph,
    partitions.total_length(),
    &partitions,
    &indels,
    settings.opt_method,
    &mut state.branch_lengths,
  )?;
  state.convergence.record(&unconverged, graph, &state.names);
  let zero_optimal_edges = if topology_ops.collapse_short_branches {
    find_zero_optimal_internal_edges(graph, &indels, &state.branch_lengths)
  } else {
    vec![]
  };
  apply_damping(
    &mut state.branch_lengths,
    &old_branch_lengths,
    settings.damping,
    iteration,
  );

  let cleanup = prune_and_merge_in_loop(
    graph,
    partitions,
    &zero_optimal_edges,
    topology_ops,
    &mut state.branch_lengths,
    &mut state.names,
    &mut state.ledger,
  )?;
  if cleanup.topology_changed {
    state.history.topology_changed();
  }
  let partitions = update_partition_rates(
    graph,
    cleanup.partitions,
    settings.partition_rates,
    &mut state.branch_lengths,
    Some(&mut state.indels),
    log,
  )?;
  Ok((partitions, cleanup.topology_changed))
}

struct LoopState {
  branch_lengths: EdgeMap<Option<f64>>,
  names: NodeNames,
  convergence: ConvergenceSummary,
  indels: EdgeIndels,
  ledger: MergeLedger,
  lh_history: Vec<LogLh>,
  history: TopologyHistory,
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
  branch_lengths: &EdgeMap<Option<f64>>,
  partitions: PartitionSet,
  indels: &EdgeIndels,
) -> Result<OptimizeIteration, Report> {
  let (partitions, marginal_lh) = partitions.marginal_update(graph, &branch_lengths_or_zero(branch_lengths))?;
  let (sparse_lh, dense_lh) = match partitions {
    PartitionSet::Sparse(_) => (marginal_lh, LogLh::ZERO),
    PartitionSet::Dense(_) => (LogLh::ZERO, marginal_lh),
  };
  let indel_lh = indels
    .recount(graph, &partitions)?
    .total_log_lh(graph, branch_lengths)?;
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
  indels: &EdgeIndels,
  branch_lengths: &EdgeMap<Option<f64>>,
) -> Vec<GraphEdgeKey> {
  graph
    .get_edges()
    .filter_map(|edge| {
      let edge_key = edge.key();
      let bl = branch_lengths[&edge_key].unwrap_or(f64::NAN);
      let is_zero_optimal = bl == 0.0 && !graph.is_leaf(edge.target());
      (is_zero_optimal && !indels.term(edge_key).requires_positive_length()).then_some(edge_key)
    })
    .collect_vec()
}

pub(crate) fn prune_and_merge_in_loop(
  graph: &mut Graph,
  mut partitions: PartitionSet,
  zero_optimal_edges: &[GraphEdgeKey],
  topology_ops: TopologyOps,
  branch_lengths: &mut EdgeMap<Option<f64>>,
  names: &mut NodeNames,
  ledger: &mut MergeLedger,
) -> Result<TopologyCleanup, Report> {
  let sparse = partitions.sparse_mut();

  let mut topology_changed = false;

  if !zero_optimal_edges.is_empty() {
    for &edge_key in zero_optimal_edges {
      if branch_lengths.contains_key(edge_key) {
        branch_lengths.insert(edge_key, Some(0.0));
      }
    }

    let groups = sparse_alignment_groups(sparse)?;
    let mut collapsed = 0_usize;
    for &edge_key in zero_optimal_edges {
      if graph.get_edge(edge_key).is_none() {
        continue;
      }
      ledger.reject_collapsed_helper(graph, names, edge_key)?;
      collapse_edge(graph, sparse, &groups, edge_key, branch_lengths)?;
      collapsed += 1;
    }

    if collapsed > 0 {
      debug!("Collapsed {collapsed} zero-optimal internal edges");
      topology_changed = true;
    }
  }

  if resolve_polytomies(graph, names, sparse, topology_ops, branch_lengths, ledger)? > 0 {
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
  branch_lengths: &mut EdgeMap<Option<f64>>,
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
