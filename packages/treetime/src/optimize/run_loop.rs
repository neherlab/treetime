use crate::ancestral::marginal::update_marginal;
use crate::optimize::args::BranchOptMethod;
use crate::optimize::iteration::{
  apply_damping, find_zero_optimal_internal_edges, prune_and_merge_in_loop, restore_branch_lengths,
  save_branch_lengths,
};
use crate::optimize::optimize_indel::{estimate_indel_rate, total_indel_log_lh};
use crate::optimize::optimize_unified::run_optimize_mixed_inner;
use crate::representation::partition::marginal_dense::PartitionMarginalDense;
use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
use crate::representation::partition::traits::{PartitionOptimizeOps, PartitionOptimizeVec};
use crate::representation::payload::ancestral::GraphAncestral;
use eyre::Report;
use itertools::{Itertools, chain};
use log::debug;
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_utils::fmt::float::float_to_significant_digits;

/// Iterative branch-length optimization with marginal reconstruction and topology cleanup.
///
/// This is the core of `run_optimize`: the I/O-free loop that runs after inputs are loaded
/// and partitions are initialized. `run_optimize` sets up partitions and calls this; tests
/// that exercise the loop should also call this directly rather than reimplementing the body.
///
/// Pre-conditions (enforced by callers, not re-checked here):
///
/// - `graph` has initial branch lengths (e.g. from [`apply_initial_guess_mode`]).
/// - Sparse partitions have compressed sequences ([`compress_sequences`]).
/// - Dense partitions have populated profiles ([`initialize_marginal`] + [`update_marginal`]).
/// - `mixed_partitions` contains both partition families (see [`collect_optimize_partitions`]).
/// - Each partition's `gtr` field is the final model (resolved before partition construction).
///
/// The indel rate is estimated once before the first iteration and held fixed
/// throughout the loop, removing the feedback path where branch-length changes
/// shift the rate estimate. When `no_indels` is true, the rate is zero and
/// the Poisson indel term drops out of both the likelihood evaluation and the
/// per-edge optimizer.
///
/// Per-iteration sequence:
///
/// 1. Run `update_marginal` on sparse and dense partitions and sum the joint
///    substitution + indel log-likelihood using the pre-computed indel rate.
/// 2. Check three stopping conditions (converged, oscillating, worsened).
/// 3. Save current branch lengths ([`save_branch_lengths`]).
/// 4. Per-edge branch-length update with the pre-computed indel rate
///    ([`run_optimize_mixed_inner`]).
/// 5. Identify internal edges the optimizer drove to exactly zero, BEFORE damping
///    blends those zeros with the old branch lengths ([`find_zero_optimal_internal_edges`]).
/// 6. Apply damping ([`apply_damping`]): convex combination of new and saved branch lengths.
/// 7. Collapse zero-optimal internal edges and merge shared mutations in resulting
///    polytomies ([`prune_and_merge_in_loop`]).
///
/// The returned [`OptimizeLoopResult`] surfaces the per-iteration likelihood history and
/// the convergence point, which are useful for tests and diagnostics but unused by the
/// production `run_optimize` wrapper.
pub fn run_optimize_loop(
  graph: &mut GraphAncestral,
  sparse_partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
  dense_partitions: &[Arc<RwLock<PartitionMarginalDense>>],
  mixed_partitions: &PartitionOptimizeVec,
  max_iter: usize,
  dp: f64,
  damping: f64,
  opt_method: BranchOptMethod,
  no_indels: bool,
) -> Result<OptimizeLoopResult, Report> {
  let indel_rate = if no_indels {
    0.0
  } else {
    estimate_indel_rate(graph, mixed_partitions)
  };

  let mut lh_history: Vec<f64> = Vec::with_capacity(max_iter);
  let mut stopped_at: Option<(usize, ConvergenceReason)> = None;
  let mut lh_prev = f64::NEG_INFINITY;
  let mut lh_prev_prev = f64::NEG_INFINITY;

  let mut best_lh = f64::NEG_INFINITY;
  let mut best_branch_lengths: Vec<f64> = Vec::new();

  for i in 0..max_iter {
    let iteration_lh = compute_iteration_likelihood(
      graph,
      sparse_partitions,
      dense_partitions,
      mixed_partitions,
      indel_rate,
      no_indels,
    )?;
    lh_history.push(iteration_lh.total_lh);

    debug!(
      "Iteration {}: likelihood {} (sparse: {}, dense: {}, indel: {}, rate: {})",
      i + 1,
      float_to_significant_digits(iteration_lh.total_lh, 7),
      float_to_significant_digits(iteration_lh.sparse_lh, 7),
      float_to_significant_digits(iteration_lh.dense_lh, 7),
      float_to_significant_digits(iteration_lh.indel_lh, 7),
      float_to_significant_digits(iteration_lh.indel_rate, 7)
    );

    if !iteration_lh.total_lh.is_finite() {
      if best_branch_lengths.len() == graph.get_edges().len() {
        restore_branch_lengths(graph, &best_branch_lengths);
        update_marginal(graph, sparse_partitions)?;
        update_marginal(graph, dense_partitions)?;
      }
      stopped_at = Some((i, ConvergenceReason::NumericalFailure));
      break;
    }

    if iteration_lh.total_lh > best_lh {
      best_lh = iteration_lh.total_lh;
      best_branch_lengths = save_branch_lengths(graph);
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
      restore_branch_lengths(graph, &best_branch_lengths);
      update_marginal(graph, sparse_partitions)?;
      update_marginal(graph, dense_partitions)?;
      stopped_at = Some((i, ConvergenceReason::Worsened));
      break;
    }

    let old_branch_lengths = save_branch_lengths(graph);
    run_optimize_mixed_inner(graph, mixed_partitions, opt_method, indel_rate, no_indels)?;

    let zero_optimal_edges = find_zero_optimal_internal_edges(graph, sparse_partitions);

    apply_damping(graph, &old_branch_lengths, damping, i);

    let topology_changed =
      prune_and_merge_in_loop(graph, sparse_partitions, dense_partitions, &zero_optimal_edges)?;
    if topology_changed {
      best_lh = f64::NEG_INFINITY;
    }

    lh_prev_prev = lh_prev;
    lh_prev = iteration_lh.total_lh;
  }

  Ok(OptimizeLoopResult { lh_history, stopped_at })
}

/// Why the optimization loop stopped early (before exhausting `max_iter`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ConvergenceReason {
  /// Successive likelihoods differ by less than `dp`: $|LH_i - LH_{i-1}| < dp$.
  Converged,
  /// Two-step likelihoods differ by less than `dp`: $|LH_i - LH_{i-2}| < dp$.
  /// Detects 2-cycles where consecutive differences exceed `dp` but the
  /// period-2 amplitude is small.
  Oscillating,
  /// Likelihood decreased from the best observed value. Branch lengths have
  /// been reverted to the best-observed state before returning.
  Worsened,
  /// Log-likelihood became NaN or infinite, indicating numerical instability
  /// in the marginal reconstruction. Branch lengths reflect the last finite state.
  NumericalFailure,
}

/// Diagnostics from [`run_optimize_loop`].
///
/// Exposed primarily for tests; `run_optimize` discards it.
#[derive(Clone, Debug, Default)]
pub struct OptimizeLoopResult {
  /// Total log-likelihood recorded at the start of each iteration, before that iteration's
  /// branch-length update. Length equals the number of iterations actually executed
  /// (including the final iteration that triggered the convergence break, if any).
  pub lh_history: Vec<f64>,

  /// Iteration index (0-based) and reason the loop stopped early.
  /// `None` if the loop exhausted `max_iter` without meeting any stopping criterion.
  pub stopped_at: Option<(usize, ConvergenceReason)>,
}

pub fn collect_optimize_partitions(
  dense_partitions: &[Arc<RwLock<PartitionMarginalDense>>],
  sparse_partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
) -> PartitionOptimizeVec {
  chain!(
    dense_partitions
      .iter()
      .cloned()
      .map(|partition| -> Arc<RwLock<dyn PartitionOptimizeOps>> { partition }),
    sparse_partitions
      .iter()
      .cloned()
      .map(|partition| -> Arc<RwLock<dyn PartitionOptimizeOps>> { partition })
  )
  .collect_vec()
}

#[derive(Clone, Copy, Debug, Default)]
struct OptimizeIterationLikelihood {
  sparse_lh: f64,
  dense_lh: f64,
  indel_lh: f64,
  total_lh: f64,
  indel_rate: f64,
}

fn compute_iteration_likelihood(
  graph: &GraphAncestral,
  sparse_partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
  dense_partitions: &[Arc<RwLock<PartitionMarginalDense>>],
  mixed_partitions: &PartitionOptimizeVec,
  indel_rate: f64,
  no_indels: bool,
) -> Result<OptimizeIterationLikelihood, Report> {
  let sparse_lh = update_marginal(graph, sparse_partitions)?;
  let dense_lh = update_marginal(graph, dense_partitions)?;
  let indel_lh = if no_indels {
    0.0
  } else {
    total_indel_log_lh(graph, mixed_partitions, indel_rate)
  };
  let total_lh = sparse_lh + dense_lh + indel_lh;

  Ok(OptimizeIterationLikelihood {
    sparse_lh,
    dense_lh,
    indel_lh,
    total_lh,
    indel_rate,
  })
}
