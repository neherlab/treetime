use crate::ancestral::marginal::update_marginal;
use crate::optimize::dispatch::initial_guess_mixed;
use crate::optimize::dispatch::run_optimize_mixed_inner;
use crate::optimize::indel::{estimate_indel_rate, total_indel_log_lh};
use crate::optimize::iteration::{apply_damping, restore_branch_lengths, save_branch_lengths};
use crate::optimize::params::{BranchOptMethod, InitialGuessMode};
use crate::optimize::topology::collapse::collapse_edge;
use crate::optimize::topology::merge_shared_mutations::merge_shared_mutation_branches;
use crate::partition::marginal_dense::PartitionMarginalDense;
use crate::partition::marginal_sparse::PartitionMarginalSparse;
use crate::partition::traits::{HasGtr, PartitionOptimizeOps, PartitionOptimizeVec};
use crate::payload::ancestral::GraphAncestral;
use eyre::Report;
use itertools::{Itertools, chain};
use log::debug;
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_graph::edge::{GraphEdge, GraphEdgeKey, HasBranchLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNode;
use treetime_utils::fmt::float::float_to_significant_digits;
use treetime_utils::make_error;

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

    let topology_changed = prune_and_merge_in_loop(graph, sparse_partitions, dense_partitions, &zero_optimal_edges)?;
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

/// Collect internal edge keys whose branch length is exactly zero after optimization,
/// provided the edge carries no substitutions or indels on any sparse partition.
///
/// `run_optimize_mixed()` sets a branch length to exactly zero through two paths:
/// the pre-dispatch `is_zero_branch_optimal()` derivative shortcut (for models
/// with proven unimodal branch-length likelihood) and the post-dispatch
/// grid-verified boundary check (for all other models, delegating to
/// `grid_search_inner` when the optimizer's point is worse than zero). This
/// function collects internal edges that reach exactly zero through either
/// path. Call this BEFORE `apply_damping()`, which would blend the zero values
/// with old branch lengths and obscure the optimizer's decision.
///
/// Only internal (non-leaf-targeting) edges are candidates: collapsing a leaf edge
/// would remove the leaf from the tree.
///
/// Leaf edges with optimizer output `0.0` are intentionally not collected here.
/// This matches v0 (`packages/legacy/treetime/treetime/treeanc.py`
/// `prune_short_branches`, which explicitly skips terminals): `apply_damping()`
/// blends the optimizer's zero with the previous positive value, producing a
/// small positive leaf branch length that converges toward zero across iterations
/// as the geometric damping weight decays. v0 also damps every non-root edge
/// (`optimize_tree_marginal`, line 1342) and never collapses leaves.
///
/// Edges that carry mutations are excluded: collapsing them would push substitutions
/// onto child edges and trigger an oscillation where merge_shared_mutation_branches
/// re-creates the node, which then optimizes to zero again.
pub fn find_zero_optimal_internal_edges(
  graph: &GraphAncestral,
  sparse_partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
) -> Vec<GraphEdgeKey> {
  graph
    .get_edges()
    .iter()
    .filter_map(|edge_ref| {
      let edge = edge_ref.read_arc();
      let bl = edge.payload().read_arc().branch_length().unwrap_or(f64::NAN);
      let target_is_leaf = graph.is_leaf(edge.target());
      if bl != 0.0 || target_is_leaf {
        return None;
      }
      let edge_key = edge.key();
      let has_mutations = sparse_partitions.iter().any(|partition| {
        let partition = partition.read_arc();
        partition
          .edges
          .get(&edge_key)
          .is_some_and(|e| !e.fitch_subs().is_empty() || !e.indels.is_empty())
      });
      (!has_mutations).then_some(edge_key)
    })
    .collect_vec()
}

/// Collapse zero-optimal internal edges and merge shared mutations in resulting polytomies.
///
/// Analogous to v0's `prune_short_branches()` inside the optimization loop:
/// edges whose optimal branch length is zero are collapsed, simplifying the tree
/// progressively across iterations. After collapsing, sibling branches in newly
/// formed polytomies that share identical substitutions are merged under new
/// internal nodes (sparse partitions only, requires discrete substitution data).
///
/// Returns true if any topology change occurred.
pub fn prune_and_merge_in_loop(
  graph: &mut GraphAncestral,
  sparse_partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
  dense_partitions: &[Arc<RwLock<PartitionMarginalDense>>],
  zero_optimal_edges: &[GraphEdgeKey],
) -> Result<bool, Report> {
  if zero_optimal_edges.is_empty() {
    return Ok(false);
  }

  // Override damped branch lengths back to zero for edges the optimizer identified
  // as zero-optimal. Damping is a convergence aid for continuous values; it should
  // not prevent collapsing degenerate edges.
  for &edge_key in zero_optimal_edges {
    if let Some(edge) = graph.get_edge(edge_key) {
      edge.write_arc().payload().write_arc().set_branch_length(Some(0.0));
    }
  }

  let mut collapsed = 0_usize;
  for &edge_key in zero_optimal_edges {
    // Edge may already be gone if a prior collapse in this batch removed it
    // (e.g., collapsing a parent also removed a child edge that was zero-optimal)
    if graph.get_edge(edge_key).is_none() {
      continue;
    }
    collapse_edge(graph, sparse_partitions, dense_partitions, edge_key)?;
    collapsed += 1;
  }

  if collapsed == 0 {
    return Ok(false);
  }

  debug!("Collapsed {collapsed} zero-optimal internal edges");

  // Merge shared mutations in newly formed polytomies (sparse only)
  if !sparse_partitions.is_empty() {
    merge_shared_mutation_branches(graph, sparse_partitions)?;
  }

  graph.build()?;
  Ok(true)
}

/// Whether a branch length is absent or NaN.
///
/// The bio crate newick parser returns `f32::NAN` for branches without
/// a `:length` suffix. After cast to f64 and wrapping in `Some`, this
/// appears as `Some(NaN)` rather than `None`. Both represent a missing
/// branch length for optimization purposes.
fn is_branch_length_missing(bl: Option<f64>) -> bool {
  bl.is_none_or(|v| v.is_nan())
}

/// Whether any edge in the graph lacks a usable branch length.
pub fn any_edge_missing_branch_length<N, E, D>(graph: &Graph<N, E, D>) -> bool
where
  N: GraphNode,
  E: GraphEdge + HasBranchLength,
  D: Send + Sync,
{
  graph
    .get_edges()
    .iter()
    .any(|edge_ref| is_branch_length_missing(edge_ref.read_arc().payload().read_arc().branch_length()))
}

/// Whether any edge that carries indels has a zero branch length.
///
/// The Poisson indel log-likelihood $\ell(t) = k \ln(\mu t) - \mu t - \ln(k!)$
/// is $-\infty$ at $t = 0$ when $k > 0$, so a zero branch length on an
/// indel-bearing edge makes the optimization objective ill-defined. In
/// addition, `estimate_indel_rate()` returns $\hat\mu = \sum k_e / \sum t_e$;
/// when every indel-bearing edge has $t_e = 0$ the numerator is positive and
/// the denominator is zero, so the estimator falls back to $\hat\mu = 0$ and
/// `poisson_indel_log_lh()` short-circuits to zero for every subsequent
/// evaluation, silently dropping the indel contribution from the likelihood
/// throughout optimization.
///
/// The `Auto` and `Always` paths bootstrap such edges to positive branch
/// lengths in `initial_guess_mixed()` before entering `run_optimize_mixed()`.
/// The `Never` path skips `initial_guess_mixed()` entirely, so it must reject
/// this configuration at validation time instead.
pub fn any_indel_edge_has_zero_branch_length<P>(graph: &GraphAncestral, partitions: &[Arc<RwLock<P>>]) -> bool
where
  P: PartitionOptimizeOps + ?Sized,
{
  graph.get_edges().iter().any(|edge_ref| {
    let edge = edge_ref.read_arc();
    let bl = edge.payload().read_arc().branch_length().unwrap_or(0.0);
    if bl != 0.0 {
      return false;
    }
    let edge_key = edge.key();
    partitions
      .iter()
      .any(|partition| partition.read_arc().edge_indel_count(edge_key) > 0)
  })
}

/// Apply the initial-guess mode to the graph.
///
/// - `Auto`: estimate only edges with missing or invalid branch lengths.
/// - `Always`: estimate all edges, overwriting existing values.
/// - `Never`: keep input branch lengths; error if any edge lacks a usable
///   value (None/NaN) or has zero branch length while carrying indels (the
///   Poisson indel log-likelihood diverges at $t = 0$ when $k > 0$).
pub fn apply_initial_guess_mode<P>(
  graph: &GraphAncestral,
  mixed_partitions: &[Arc<RwLock<P>>],
  mode: InitialGuessMode,
  no_indels: bool,
) -> Result<(), Report>
where
  P: PartitionOptimizeOps + ?Sized,
{
  match mode {
    InitialGuessMode::Auto => initial_guess_mixed(graph, mixed_partitions, false),
    InitialGuessMode::Always => initial_guess_mixed(graph, mixed_partitions, true),
    InitialGuessMode::Never => {
      if any_edge_missing_branch_length(graph) {
        return make_error!(
          "--branch-length-initial-guess=never requires all edges to have finite branch lengths, \
           but some edges have missing or NaN values. \
           Use 'auto' to fill in missing values or 'always' to recompute all branch lengths"
        );
      }
      if !no_indels && any_indel_edge_has_zero_branch_length(graph, mixed_partitions) {
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

/// Normalize substitution rates across partitions after GTR inference.
///
/// Each inferred GTR model has a rate `mu` (expected substitutions per site per branch-length
/// unit). When partitions differ in rate, we rescale so the length-weighted average mu is 1:
///
///   total_average = Σ(partition.length * partition.gtr.mu) / Σ(partition.length)
///
/// Each partition's `gtr.mu` is divided by `total_average`, and every branch length in the tree
/// is scaled by `total_average`. The value `mu * t` (expected substitutions) is preserved,
/// but now the average rate across all partitions equals 1, making branch lengths directly
/// interpretable as substitutions per site.
pub fn normalize_partition_rates<P: HasGtr>(graph: &GraphAncestral, partitions: &[Arc<RwLock<P>>]) {
  let total_length: usize = partitions.iter().map(|p| p.read_arc().sequence_length()).sum();

  if total_length == 0 {
    return;
  }

  let weighted_rate: f64 = partitions.iter().map(|p| p.read_arc().weighted_rate()).sum();

  let total_average = weighted_rate / total_length as f64;

  if total_average == 0.0 {
    return;
  }

  for partition in partitions {
    partition.write_arc().normalize_rate(total_average);
  }

  for edge_ref in graph.get_edges() {
    let mut edge = edge_ref.write_arc().payload().write_arc();
    if let Some(bl) = edge.branch_length() {
      edge.set_branch_length(Some(bl * total_average));
    }
  }
}
