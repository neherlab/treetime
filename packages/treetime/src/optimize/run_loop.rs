use crate::ancestral::marginal::branch_lengths_or_zero;
use crate::ancestral::pipeline::{DenseReconstruction, SparseReconstruction};
use crate::gtr::gtr::GTR;
use crate::optimize::branch_length::invalid_branch_length_descriptions;
use crate::optimize::dispatch::initial_guess_mixed;
use crate::optimize::dispatch::run_optimize_mixed_inner;
use crate::optimize::gather::{gather_edge_contributions, gather_edge_indel_counts, total_sequence_length};
use crate::optimize::indel::{estimate_indel_rate, total_indel_log_lh};
use crate::optimize::iteration::apply_damping;
use crate::optimize::params::ExistingBranchLengths;
use crate::optimize::params::{BranchOptMethod, InitialGuessMode, TopologyOps};
use crate::optimize::topology::collapse::collapse_edge;
use crate::optimize::topology::resolve_polytomy::resolve_polytomies;
use crate::partition::marginal::dense::partition::{DenseMarginalEdges, PartitionMarginalDense};
use crate::partition::marginal::shared::reconcile::{live_node_keys, reconcile_node_states};
use crate::partition::marginal::sparse::partition::{PartitionMarginalSparse, SparseMarginalEdges};
use crate::partition::storage::dense::DenseNodeState;
use crate::partition::storage::sparse::SparseNodeState;
use eyre::Report;
use itertools::{Itertools, izip};
use log::{debug, warn};
use std::collections::BTreeMap;
use treetime_graph::assign_node_names::assign_node_names;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::LogLh;
use treetime_utils::fmt::float::float_to_significant_digits;
use treetime_utils::make_error;

pub(crate) fn run_optimize_loop(
  graph: &mut Graph,
  sparse_partitions: Vec<SparseReconstruction>,
  dense_partitions: Vec<DenseReconstruction>,
  max_iter: usize,
  dp: f64,
  damping: f64,
  opt_method: BranchOptMethod,
  no_indels: bool,
  topology_ops: TopologyOps,
  branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Result<OptimizeLoopResult, Report> {
  let mut names = names.clone();
  let mut branch_lengths = branch_lengths;
  let mut sparse_partitions = sparse_partitions;
  let mut dense_partitions = dense_partitions;

  let indel_rate = if no_indels {
    0.0
  } else {
    let indel_counts = gather_edge_indel_counts(graph, &dense_partitions, &sparse_partitions);
    estimate_indel_rate(graph, &indel_counts, &branch_lengths)
  };

  let mut lh_history: Vec<LogLh> = Vec::with_capacity(max_iter);
  let mut stopped_at: Option<(usize, ConvergenceReason)> = None;
  let mut lh_prev = LogLh::IMPOSSIBLE;
  let mut lh_prev_prev = LogLh::IMPOSSIBLE;

  let mut best_lh = LogLh::IMPOSSIBLE;
  let mut best_branch_lengths: Option<BTreeMap<GraphEdgeKey, Option<f64>>> = None;

  for i in 0..max_iter {
    let iteration = compute_iteration(
      graph,
      &branch_lengths,
      sparse_partitions,
      dense_partitions,
      indel_rate,
      no_indels,
    )?;
    sparse_partitions = iteration.sparse_partitions;
    dense_partitions = iteration.dense_partitions;
    let iteration_lh = iteration.likelihood;
    lh_history.push(iteration_lh.total_lh);

    debug!(
      "Iteration {}: likelihood {} (sparse: {}, dense: {}, indel: {}, rate: {})",
      i + 1,
      float_to_significant_digits(iteration_lh.total_lh.value(), 7),
      float_to_significant_digits(iteration_lh.sparse_lh.value(), 7),
      float_to_significant_digits(iteration_lh.dense_lh.value(), 7),
      float_to_significant_digits(iteration_lh.indel_lh.value(), 7),
      float_to_significant_digits(iteration_lh.indel_rate, 7)
    );

    if !iteration_lh.total_lh.value().is_finite() {
      if let Some(best) = &best_branch_lengths {
        branch_lengths = best.clone();
        let marginal_bl = branch_lengths_or_zero(&branch_lengths);
        (sparse_partitions, _) = marginal_update_sparse(graph, &marginal_bl, sparse_partitions)?;
        (dense_partitions, _) = marginal_update_dense(graph, &marginal_bl, dense_partitions)?;
      }
      stopped_at = Some((i, ConvergenceReason::NumericalFailure));
      break;
    }

    if iteration_lh.total_lh > best_lh {
      best_lh = iteration_lh.total_lh;
      best_branch_lengths = Some(branch_lengths.clone());
    }

    if (iteration_lh.total_lh - lh_prev).abs() < dp.abs() {
      stopped_at = Some((i, ConvergenceReason::Converged));
      break;
    }

    if i >= 2 && (iteration_lh.total_lh - lh_prev_prev).abs() < dp.abs() {
      stopped_at = Some((i, ConvergenceReason::Oscillating));
      break;
    }

    if i >= 2 && iteration_lh.total_lh < lh_prev && lh_prev >= best_lh {
      if let Some(best) = &best_branch_lengths {
        branch_lengths = best.clone();
        let marginal_bl = branch_lengths_or_zero(&branch_lengths);
        (sparse_partitions, _) = marginal_update_sparse(graph, &marginal_bl, sparse_partitions)?;
        (dense_partitions, _) = marginal_update_dense(graph, &marginal_bl, dense_partitions)?;
      }
      stopped_at = Some((i, ConvergenceReason::Worsened));
      break;
    }

    let old_branch_lengths = branch_lengths.clone();
    {
      let total_length = total_sequence_length(&dense_partitions, &sparse_partitions);
      let contributions = gather_edge_contributions(graph, &dense_partitions, &sparse_partitions)?;
      let indel_counts = gather_edge_indel_counts(graph, &dense_partitions, &sparse_partitions);
      run_optimize_mixed_inner(
        graph,
        total_length,
        &contributions,
        &indel_counts,
        opt_method,
        indel_rate,
        no_indels,
        &mut branch_lengths,
      )?;
    }

    let zero_optimal_edges = if topology_ops.collapse_short_branches {
      find_zero_optimal_internal_edges(graph, &sparse_partitions, &branch_lengths)
    } else {
      vec![]
    };

    apply_damping(&mut branch_lengths, &old_branch_lengths, damping, i);

    let cleanup = prune_and_merge_in_loop(
      graph,
      sparse_partitions,
      dense_partitions,
      &zero_optimal_edges,
      topology_ops,
      &mut branch_lengths,
      &mut names,
    )?;
    sparse_partitions = cleanup.sparse_partitions;
    dense_partitions = cleanup.dense_partitions;
    if cleanup.topology_changed {
      best_lh = LogLh::IMPOSSIBLE;
      best_branch_lengths = None;
    }

    lh_prev_prev = lh_prev;
    lh_prev = iteration_lh.total_lh;
  }

  Ok(OptimizeLoopResult {
    sparse_partitions,
    dense_partitions,
    branch_lengths,
    names,
    lh_history,
    stopped_at,
  })
}

#[derive(Clone, Debug, Default)]
pub(crate) struct OptimizeLoopResult {
  pub sparse_partitions: Vec<SparseReconstruction>,

  pub dense_partitions: Vec<DenseReconstruction>,

  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,

  pub names: BTreeMap<GraphNodeKey, Option<String>>,

  #[allow(dead_code, reason = "read only in test code")]
  pub lh_history: Vec<LogLh>,

  #[allow(dead_code, reason = "read only in test code")]
  pub stopped_at: Option<(usize, ConvergenceReason)>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ConvergenceReason {
  Converged,
  Oscillating,
  Worsened,
  NumericalFailure,
}

fn compute_iteration(
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  sparse_partitions: Vec<SparseReconstruction>,
  dense_partitions: Vec<DenseReconstruction>,
  indel_rate: f64,
  no_indels: bool,
) -> Result<OptimizeIteration, Report> {
  let marginal_bl = branch_lengths_or_zero(branch_lengths);
  let (sparse_partitions, sparse_lh) = marginal_update_sparse(graph, &marginal_bl, sparse_partitions)?;
  let (dense_partitions, dense_lh) = marginal_update_dense(graph, &marginal_bl, dense_partitions)?;
  let indel_lh = if no_indels {
    LogLh::ZERO
  } else {
    let indel_counts = gather_edge_indel_counts(graph, &dense_partitions, &sparse_partitions);
    total_indel_log_lh(graph, &indel_counts, branch_lengths, indel_rate)?
  };
  let total_lh = sparse_lh + dense_lh + indel_lh;

  Ok(OptimizeIteration {
    sparse_partitions,
    dense_partitions,
    likelihood: OptimizeIterationLikelihood {
      sparse_lh,
      dense_lh,
      indel_lh,
      total_lh,
      indel_rate,
    },
  })
}

pub(crate) fn marginal_update_sparse(
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  sparse: Vec<SparseReconstruction>,
) -> Result<(Vec<SparseReconstruction>, LogLh), Report> {
  sparse
    .into_iter()
    .try_fold((Vec::new(), LogLh::ZERO), |(mut updated, total), family| {
      let (family, log_lh) = family.marginal_update(graph, branch_lengths)?;
      updated.push(family);
      Ok((updated, total + log_lh))
    })
}

pub(crate) fn marginal_update_dense(
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  dense: Vec<DenseReconstruction>,
) -> Result<(Vec<DenseReconstruction>, LogLh), Report> {
  dense
    .into_iter()
    .try_fold((Vec::new(), LogLh::ZERO), |(mut updated, total), family| {
      let (family, log_lh) = family.marginal_update(graph, branch_lengths)?;
      updated.push(family);
      Ok((updated, total + log_lh))
    })
}

struct OptimizeIteration {
  sparse_partitions: Vec<SparseReconstruction>,
  dense_partitions: Vec<DenseReconstruction>,
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
  sparse_partitions: &[SparseReconstruction],
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Vec<GraphEdgeKey> {
  graph
    .get_edges()
    .filter_map(|edge_ref| {
      let edge = edge_ref;
      let bl = branch_lengths[&edge.key()].unwrap_or(f64::NAN);
      let target_is_leaf = graph.is_leaf(edge.target());
      if bl != 0.0 || target_is_leaf {
        return None;
      }
      let edge_key = edge.key();
      let has_mutations = sparse_partitions.iter().any(|family| {
        family
          .partition
          .obs_edges
          .get(&edge_key)
          .is_some_and(|e| !e.fitch_subs().is_empty() || !e.indels.is_empty())
      });
      (!has_mutations).then_some(edge_key)
    })
    .collect_vec()
}

pub(crate) fn prune_and_merge_in_loop(
  graph: &mut Graph,
  sparse_partitions: Vec<SparseReconstruction>,
  dense_partitions: Vec<DenseReconstruction>,
  zero_optimal_edges: &[GraphEdgeKey],
  topology_ops: TopologyOps,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  names: &mut BTreeMap<GraphNodeKey, Option<String>>,
) -> Result<TopologyCleanup, Report> {
  let (mut sparse_obs, sparse_gtrs, mut sparse_node_states, sparse_edges): (Vec<_>, Vec<_>, Vec<_>, Vec<_>) =
    sparse_partitions
      .into_iter()
      .map(|family| (family.partition, family.gtr, family.node_states, family.edges))
      .multiunzip();
  let (dense_obs, dense_gtrs, dense_node_states, dense_edges): (Vec<_>, Vec<_>, Vec<_>, Vec<_>) = dense_partitions
    .into_iter()
    .map(|family| (family.partition, family.gtr, family.node_states, family.edges))
    .multiunzip();

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
      collapse_edge(graph, &mut sparse_obs, edge_key, branch_lengths)?;
      collapsed += 1;
    }

    if collapsed > 0 {
      debug!("Collapsed {collapsed} zero-optimal internal edges");
      topology_changed = true;
    }
  }

  if resolve_polytomies(
    graph,
    &mut sparse_obs,
    &mut sparse_node_states,
    topology_ops,
    branch_lengths,
  )? > 0
  {
    topology_changed = true;
  }

  if !topology_changed {
    return Ok(TopologyCleanup {
      sparse_partitions: join_sparse(sparse_obs, sparse_gtrs, sparse_node_states, sparse_edges),
      dense_partitions: join_dense(dense_obs, dense_gtrs, dense_node_states, dense_edges),
      topology_changed,
    });
  }

  graph.build()?;
  *names = assign_node_names(std::mem::take(names), graph)?;
  let live_nodes = live_node_keys(graph);
  let sparse_partitions = izip!(sparse_obs, sparse_gtrs, sparse_node_states)
    .map(|(partition, gtr, node_states)| {
      SparseReconstruction::seeded(
        partition,
        gtr,
        reconcile_node_states(node_states, &live_nodes, SparseNodeState::empty),
      )
    })
    .collect_vec();
  let dense_partitions = izip!(dense_obs, dense_gtrs, dense_node_states)
    .map(|(partition, gtr, node_states)| {
      DenseReconstruction::seeded(
        partition,
        gtr,
        reconcile_node_states(node_states, &live_nodes, DenseNodeState::empty),
      )
    })
    .collect_vec();

  Ok(TopologyCleanup {
    sparse_partitions,
    dense_partitions,
    topology_changed,
  })
}

pub(crate) struct TopologyCleanup {
  pub sparse_partitions: Vec<SparseReconstruction>,
  pub dense_partitions: Vec<DenseReconstruction>,
  pub topology_changed: bool,
}

fn join_sparse(
  partitions: Vec<PartitionMarginalSparse>,
  gtrs: Vec<GTR>,
  node_states: Vec<BTreeMap<GraphNodeKey, SparseNodeState>>,
  edges: Vec<SparseMarginalEdges>,
) -> Vec<SparseReconstruction> {
  izip!(partitions, gtrs, node_states, edges)
    .map(|(partition, gtr, node_states, edges)| SparseReconstruction {
      partition,
      gtr,
      node_states,
      edges,
    })
    .collect_vec()
}

fn join_dense(
  partitions: Vec<PartitionMarginalDense>,
  gtrs: Vec<GTR>,
  node_states: Vec<BTreeMap<GraphNodeKey, DenseNodeState>>,
  edges: Vec<DenseMarginalEdges>,
) -> Vec<DenseReconstruction> {
  izip!(partitions, gtrs, node_states, edges)
    .map(|(partition, gtr, node_states, edges)| DenseReconstruction {
      partition,
      gtr,
      node_states,
      edges,
    })
    .collect_vec()
}

#[expect(
  clippy::too_many_arguments,
  reason = "each argument is an independent input of this step; a parameter struct would be built only for this call"
)]
pub(crate) fn apply_initial_guess_mode(
  graph: &Graph,
  total_length: usize,
  indel_counts: &BTreeMap<GraphEdgeKey, usize>,
  sub_counts: &BTreeMap<GraphEdgeKey, usize>,
  effective_lengths: &BTreeMap<GraphEdgeKey, usize>,
  mode: InitialGuessMode,
  no_indels: bool,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Result<(), Report> {
  let invalid_branch_lengths = invalid_branch_length_descriptions(graph, branch_lengths, names)?;
  if let Some(message) = invalid_branch_length_warning(&invalid_branch_lengths) {
    warn!("{message}");
  }

  match mode {
    InitialGuessMode::Auto => initial_guess_mixed(
      graph,
      total_length,
      indel_counts,
      sub_counts,
      effective_lengths,
      ExistingBranchLengths::Keep,
      no_indels,
      branch_lengths,
    ),
    InitialGuessMode::Always => initial_guess_mixed(
      graph,
      total_length,
      indel_counts,
      sub_counts,
      effective_lengths,
      ExistingBranchLengths::Overwrite,
      no_indels,
      branch_lengths,
    ),
    InitialGuessMode::Never => {
      if !invalid_branch_lengths.is_empty() {
        return make_error!(
          "--branch-length-initial-guess=never requires every edge to have a finite, non-negative branch length. \
           Invalid edges:\n  {}\n\
           Use 'auto' to fill in missing values or 'always' to recompute all branch lengths",
          invalid_branch_lengths.join("\n  ")
        );
      }
      if !no_indels && any_indel_edge_has_zero_branch_length(graph, indel_counts, branch_lengths) {
        return make_error!(
          "--branch-length-initial-guess=never requires non-zero branch lengths on edges that carry indels, \
           but some indel-bearing edges have branch length zero. \
           The Poisson indel log-likelihood diverges at t=0 when k>0, so zero is not a valid input there. \
           Use 'auto' or 'always' to bootstrap these edges, or provide positive branch lengths"
        );
      }
      Ok(())
    },
  }
}

pub(crate) fn any_indel_edge_has_zero_branch_length(
  graph: &Graph,
  indel_counts: &BTreeMap<GraphEdgeKey, usize>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
) -> bool {
  graph.get_edges().any(|edge_ref| {
    let edge = edge_ref;
    let edge_key = edge.key();
    let bl = branch_lengths[&edge_key].unwrap_or(0.0);
    if bl != 0.0 {
      return false;
    }
    indel_counts[&edge_key] > 0
  })
}

pub(super) fn invalid_branch_length_warning(invalid_branch_lengths: &[String]) -> Option<String> {
  (!invalid_branch_lengths.is_empty()).then(|| {
    format!(
      "Input tree contains invalid branch lengths (expected finite values >= 0):\n  {}\n\
       Use --branch-length-initial-guess=auto to replace invalid values or \
       --branch-length-initial-guess=always to recompute all branch lengths",
      invalid_branch_lengths.join("\n  ")
    )
  })
}

#[allow(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
pub(crate) fn normalize_partition_rates(
  partitions: &mut [(usize, &mut GTR)],
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
) {
  let total_length: usize = partitions.iter().map(|(length, _)| *length).sum();

  if total_length == 0 {
    return;
  }

  let weighted_rate: f64 = partitions.iter().map(|(length, gtr)| *length as f64 * gtr.mu).sum();

  let total_average = weighted_rate / total_length as f64;

  if total_average == 0.0 {
    return;
  }

  for (_, gtr) in partitions.iter_mut() {
    gtr.mu /= total_average;
  }

  for value in branch_lengths.values_mut().flatten() {
    *value *= total_average;
  }
}
