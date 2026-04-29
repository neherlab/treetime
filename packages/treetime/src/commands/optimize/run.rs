use crate::alphabet::alphabet::Alphabet;
use crate::commands::ancestral::fitch::{compress_sequences, get_common_length};
use crate::commands::ancestral::marginal::{initialize_marginal, update_marginal};
use crate::commands::optimize::args::{BranchOptMethod, InitialGuessMode, TreetimeOptimizeArgs};
use crate::commands::optimize::optimize_indel::{estimate_indel_rate, total_indel_log_lh};
use crate::commands::optimize::optimize_unified::{initial_guess_mixed, run_optimize_mixed_inner};
use crate::commands::optimize::partition_ops::{PartitionOptimizeOps, PartitionOptimizeVec};
use crate::commands::prune::run::merge_shared_mutation_branches;
use crate::gtr::get_gtr::{GtrModelName, JC69Params, get_gtr_dense, get_gtr_sparse, jc69, write_gtr_json};
use crate::representation::algo::infer_dense::infer_dense;
use crate::representation::algo::topology_cleanup::collapse::collapse_edge;
use crate::representation::partition::marginal_dense::PartitionMarginalDense;
use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
use crate::representation::partition::traits::PartitionBranchOps;
use crate::representation::payload::ancestral::{GraphAncestral, annotate_branch_mutations};
use eyre::Report;
use itertools::{Itertools, chain, izip};
use log::debug;
use maplit::btreemap;
use num_traits::pow::pow;
use parking_lot::RwLock;
use serde::Serialize;
use std::path::Path;
use std::sync::Arc;
use treetime_graph::edge::{GraphEdge, GraphEdgeKey, HasBranchLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNode;
use treetime_io::fasta::read_many_fasta;
use treetime_io::nex::{NexWriteOptions, nex_write_file};
use treetime_io::nwk::{EdgeToNwk, NodeToNwk, NwkWriteOptions, nwk_read_file, nwk_write_file};
use treetime_primitives::seq;
use treetime_utils::fmt::float::float_to_significant_digits;
use treetime_utils::make_error;

#[derive(Clone, Debug, Default)]
pub struct TreetimeOptimizeParams {
  pub sample_from_profile: bool,
  pub fixed_pi: bool,
}

/// Normalize substitution rates across partitions after GTR inference.
///
/// Each inferred GTR model has a rate `mu` (expected substitutions per site per branch-length
/// unit). When partitions differ in rate, we rescale so the length-weighted average mu is 1:
///
///   total_average = Σ(partition.length × partition.gtr.mu) / Σ(partition.length)
///
/// Each partition's `gtr.mu` is divided by `total_average`, and every branch length in the tree
/// is multiplied by `total_average`. The product `mu × t` (expected substitutions) is preserved,
/// but now the average rate across all partitions equals 1, making branch lengths directly
/// interpretable as substitutions per site.
fn normalize_partition_rates(
  graph: &GraphAncestral,
  sparse_partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
  dense_partitions: &[Arc<RwLock<PartitionMarginalDense>>],
) {
  let total_length: usize = sparse_partitions.iter().map(|p| p.read_arc().length).sum::<usize>()
    + dense_partitions.iter().map(|p| p.read_arc().length).sum::<usize>();

  if total_length == 0 {
    return;
  }

  let weighted_rate: f64 = sparse_partitions
    .iter()
    .map(|p| {
      let p = p.read_arc();
      p.length as f64 * p.gtr.mu
    })
    .sum::<f64>()
    + dense_partitions
      .iter()
      .map(|p| {
        let p = p.read_arc();
        p.length as f64 * p.gtr.mu
      })
      .sum::<f64>();

  let total_average = weighted_rate / total_length as f64;

  if total_average == 0.0 {
    return;
  }

  for partition in sparse_partitions {
    partition.write_arc().gtr.mu /= total_average;
  }
  for partition in dense_partitions {
    partition.write_arc().gtr.mu /= total_average;
  }

  for edge_ref in graph.get_edges() {
    let mut edge = edge_ref.write_arc().payload().write_arc();
    if let Some(bl) = edge.branch_length() {
      edge.set_branch_length(Some(bl * total_average));
    }
  }
}

pub fn run_optimize(args: &TreetimeOptimizeArgs) -> Result<(), Report> {
  let TreetimeOptimizeArgs {
    input_fastas,
    tree,
    alphabet,
    model_name,
    dense,
    outdir,
    max_iter,
    dp,
    damping,
    branch_length_initial_guess,
    opt_method,
    no_indels,
  } = args;

  if !(0.0..1.0).contains(damping) {
    return make_error!("--damping must be in [0.0, 1.0), got {damping}");
  }

  let dense = dense.unwrap_or_else(infer_dense);
  let alphabet = Alphabet::new(alphabet.unwrap_or_default())?;
  let aln = read_many_fasta(input_fastas, &alphabet)?;
  let mut graph: GraphAncestral = nwk_read_file(tree)?;

  let sparse_partitions: Vec<Arc<RwLock<PartitionMarginalSparse>>> = if !dense {
    #[allow(clippy::iter_on_single_items)]
    let partitions = [PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?, // FIXME: dummy temporary gtr should not be needed here
      alphabet: alphabet.clone(),
      length: get_common_length(&aln)?,
      root_sequence: seq![],
      nodes: btreemap! {},
      edges: btreemap! {},
    }]
    .into_iter()
    .map(|p| Arc::new(RwLock::new(p)))
    .collect_vec();

    compress_sequences(&graph, &partitions, &aln)?;

    // FIXME: chicken & egg problem: to get a gtr we need partitions, to get partitions we need a gtr
    // FIXME: spaghetti code: dummy gtr is replaced by real gtr here
    for partition in &partitions {
      let gtr = get_gtr_sparse(model_name, partition, &graph)?;
      partition.write_arc().gtr = gtr;
    }

    partitions
  } else {
    vec![]
  };

  let dense_partitions: Vec<Arc<RwLock<PartitionMarginalDense>>> = if dense {
    #[allow(clippy::iter_on_single_items)]
    let partitions = [PartitionMarginalDense {
      index: 0,
      gtr: jc69(JC69Params::default())?, // FIXME: dummy temporary gtr should not be needed here
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }]
    .into_iter()
    .map(|p| Arc::new(RwLock::new(p)))
    .collect_vec();

    partitions
  } else {
    vec![]
  };

  // Initialize marginal data for whichever partitions are active
  update_marginal(&graph, &sparse_partitions)?;
  if !dense_partitions.is_empty() {
    initialize_marginal(&graph, &dense_partitions, &aln)?;

    // FIXME: chicken & egg problem: to get a gtr we need partitions, to get partitions we need a gtr
    // Dense GTR requires profiles populated via initialize_marginal + update_marginal.
    // Run first marginal pass with dummy GTR, then replace with real model.
    update_marginal(&graph, &dense_partitions)?;
    for partition in &dense_partitions {
      let gtr = get_gtr_dense(model_name, partition, &graph)?;
      partition.write_arc().gtr = gtr;
    }

    // Re-run with inferred GTR so node posteriors reflect the real model.
    update_marginal(&graph, &dense_partitions)?;
  }

  let mixed_partitions = collect_optimize_partitions(&dense_partitions, &sparse_partitions);

  // When GTR is inferred, normalize substitution rates so the length-weighted average mu equals 1
  // and rescale branch lengths accordingly. This makes branch lengths interpretable as expected
  // substitutions per site across all partitions, while preserving relative rates between partitions.
  if *model_name == GtrModelName::Infer {
    normalize_partition_rates(&graph, &sparse_partitions, &dense_partitions);
  }

  for partition in &sparse_partitions {
    write_gtr_json(&partition.read_arc().gtr, *model_name, outdir, None)?;
  }
  for partition in &dense_partitions {
    write_gtr_json(&partition.read_arc().gtr, *model_name, outdir, None)?;
  }

  apply_initial_guess_mode(&graph, &mixed_partitions, *branch_length_initial_guess, *no_indels)?;

  run_optimize_loop(
    &mut graph,
    &sparse_partitions,
    &dense_partitions,
    &mixed_partitions,
    *max_iter,
    *dp,
    *damping,
    *opt_method,
    *no_indels,
  )?;

  // Re-run marginal to populate subs_ml after the loop (the last iteration
  // may have changed topology, clearing ML subs on affected edges).
  update_marginal(&graph, &sparse_partitions)?;
  update_marginal(&graph, &dense_partitions)?;

  // Annotate mutations from whichever partition type is active
  let branch_ops: Vec<Arc<RwLock<dyn PartitionBranchOps>>> = if dense {
    dense_partitions
      .iter()
      .cloned()
      .map(|p| -> Arc<RwLock<dyn PartitionBranchOps>> { p })
      .collect_vec()
  } else {
    sparse_partitions
      .iter()
      .cloned()
      .map(|p| -> Arc<RwLock<dyn PartitionBranchOps>> { p })
      .collect_vec()
  };
  annotate_branch_mutations(&graph, &branch_ops)?;

  write_graph(outdir, &graph)?;
  Ok(())
}

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
/// - Each partition's `gtr` field is the final model (dummy GTRs have been replaced).
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
  // Compute the indel rate once before entering the loop, rather than
  // re-estimating each iteration. This removes the feedback path where
  // branch-length changes shift the rate estimate, which shifts
  // branch-length targets. When --no-indels is set, the rate is zero
  // and the Poisson indel term drops out entirely.
  let indel_rate = if no_indels {
    0.0
  } else {
    estimate_indel_rate(graph, mixed_partitions)
  };

  let mut lh_history: Vec<f64> = Vec::with_capacity(max_iter);
  let mut stopped_at: Option<(usize, ConvergenceReason)> = None;
  let mut lh_prev = f64::NEG_INFINITY;
  let mut lh_prev_prev = f64::NEG_INFINITY;

  // Best-observed state for revert on the worsened condition.
  let mut best_lh = f64::NEG_INFINITY;
  let mut best_branch_lengths: Vec<f64> = Vec::new();

  for i in 0..max_iter {
    let iteration_lh = compute_iteration_likelihood(graph, sparse_partitions, dense_partitions, mixed_partitions, indel_rate, no_indels)?;
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

    // NaN/Inf from numerical instability silently bypasses all convergence
    // checks (NaN < x is false for all x under IEEE 754). Record the failure
    // so callers can distinguish numerical breakdown from normal termination.
    if !iteration_lh.total_lh.is_finite() {
      if best_branch_lengths.len() == graph.get_edges().len() {
        restore_branch_lengths(graph, &best_branch_lengths);
        update_marginal(graph, sparse_partitions)?;
        update_marginal(graph, dense_partitions)?;
      }
      stopped_at = Some((i, ConvergenceReason::NumericalFailure));
      break;
    }

    // Track best-observed state for potential revert.
    if iteration_lh.total_lh > best_lh {
      best_lh = iteration_lh.total_lh;
      best_branch_lengths = save_branch_lengths(graph);
    }

    // Three orthogonal stopping conditions, checked in order of specificity.
    //
    // 1. Converged: successive likelihoods within dp (standard EM criterion).
    if (iteration_lh.total_lh - lh_prev).abs() < dp.abs() {
      stopped_at = Some((i, ConvergenceReason::Converged));
      break;
    }

    // 2. Oscillating: two-step likelihoods within dp (detects 2-cycles where
    //    consecutive differences exceed dp but the period-2 amplitude is small).
    //    Requires at least 2 prior measurements (i >= 2).
    if i >= 2 && (iteration_lh.total_lh - lh_prev_prev).abs() < dp.abs() {
      stopped_at = Some((i, ConvergenceReason::Oscillating));
      break;
    }

    // 3. Worsened: likelihood decreased from the best-observed value.
    //    Revert branch lengths to the best state and re-run marginal
    //    reconstruction so partition posteriors are consistent with the
    //    restored branch lengths before returning.
    //    The `lh_prev >= best_lh` guard prevents firing during normal
    //    non-monotone transients where the overall trend is upward.
    //    Requires i >= 2 to allow the initial transient.
    if i >= 2 && iteration_lh.total_lh < lh_prev && lh_prev >= best_lh {
      restore_branch_lengths(graph, &best_branch_lengths);
      update_marginal(graph, sparse_partitions)?;
      update_marginal(graph, dense_partitions)?;
      stopped_at = Some((i, ConvergenceReason::Worsened));
      break;
    }

    let old_branch_lengths = save_branch_lengths(graph);
    run_optimize_mixed_inner(graph, mixed_partitions, opt_method, indel_rate, no_indels)?;

    // Identify zero-optimal internal edges BEFORE damping blends them with old values
    let zero_optimal_edges = find_zero_optimal_internal_edges(graph, sparse_partitions);

    apply_damping(graph, &old_branch_lengths, damping, i);

    // Collapse zero-optimal edges and merge shared mutations in resulting polytomies.
    // Topology changes invalidate the saved best_branch_lengths (edge count/order changed).
    let topology_changed = prune_and_merge_in_loop(graph, sparse_partitions, dense_partitions, &zero_optimal_edges)?;
    if topology_changed {
      best_lh = f64::NEG_INFINITY;
    }

    lh_prev_prev = lh_prev;
    lh_prev = iteration_lh.total_lh;
  }

  Ok(OptimizeLoopResult { lh_history, stopped_at })
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

pub(crate) fn collect_optimize_partitions(
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

/// Save current branch lengths for all edges in graph traversal order.
pub(crate) fn save_branch_lengths<N, E, D>(graph: &Graph<N, E, D>) -> Vec<f64>
where
  N: GraphNode,
  E: GraphEdge + HasBranchLength,
  D: Send + Sync,
{
  graph
    .get_edges()
    .iter()
    .map(|edge_ref| edge_ref.read_arc().payload().read_arc().branch_length().unwrap_or(0.0))
    .collect_vec()
}

/// Restore previously saved branch lengths to all edges in graph traversal order.
///
/// Inverse of [`save_branch_lengths`]. The saved vector must have been produced
/// from the same graph (same edge count and order). Used by the worsened-condition
/// revert in [`run_optimize_loop`].
pub(crate) fn restore_branch_lengths<N, E, D>(graph: &Graph<N, E, D>, saved: &[f64])
where
  N: GraphNode,
  E: GraphEdge + HasBranchLength,
  D: Send + Sync,
{
  debug_assert_eq!(
    graph.get_edges().len(),
    saved.len(),
    "restore_branch_lengths: edge count changed between save and restore"
  );
  for (edge_ref, &bl) in izip!(graph.get_edges(), saved) {
    edge_ref.write_arc().payload().write_arc().set_branch_length(Some(bl));
  }
}

/// Minimum fraction of the old branch length retained at any iteration.
///
/// Without a floor, exponential damping $d^{i+1}$ decays to effectively zero
/// at high iteration counts (e.g. $0.75^{20} \approx 0.003$). On datasets where
/// the sparse variable/fixed position boundary oscillates, fully undamped late
/// iterations amplify the discrete jump. The floor ensures at least 1% of the
/// old value is retained, bridging the discontinuity at all iteration counts.
pub(crate) const DAMPING_FLOOR: f64 = 0.01;

/// Blend optimized branch lengths with saved old values using exponential damping.
///
/// At iteration `i` (0-based), each branch length becomes:
///   bl = bl_optimized * (1 - damping_factor) + bl_old * damping_factor
///
/// where `damping_factor = max(damping^(i+1), DAMPING_FLOOR)`.
///
/// When `damping == 0.0`, damping_factor = 0 and the optimized value is kept unchanged.
/// Early iterations take conservative steps; later iterations approach the full update
/// but never go below the `DAMPING_FLOOR` weight on the old value.
pub(crate) fn apply_damping<N, E, D>(graph: &Graph<N, E, D>, old_branch_lengths: &[f64], damping: f64, iteration: usize)
where
  N: GraphNode,
  E: GraphEdge + HasBranchLength,
  D: Send + Sync,
{
  if damping == 0.0 {
    return;
  }
  debug_assert_eq!(
    graph.get_edges().len(),
    old_branch_lengths.len(),
    "apply_damping: edge count changed between save and apply"
  );
  let damping_factor = pow(damping, iteration + 1).max(DAMPING_FLOOR);
  let new_weight = 1.0 - damping_factor;
  for (edge_ref, &old_bl) in izip!(graph.get_edges(), old_branch_lengths) {
    let mut edge = edge_ref.write_arc().payload().write_arc();
    let optimized_bl = edge.branch_length().unwrap_or(0.0);
    let damped_bl = optimized_bl * new_weight + old_bl * damping_factor;
    edge.set_branch_length(Some(damped_bl));
  }
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
/// onto child edges and potentially trigger an oscillation where merge_shared_mutation_branches
/// re-creates the node, which then optimizes to zero again.
pub(crate) fn find_zero_optimal_internal_edges(
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
pub(crate) fn prune_and_merge_in_loop(
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
pub(crate) fn any_edge_missing_branch_length<N, E, D>(graph: &Graph<N, E, D>) -> bool
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
pub(crate) fn any_indel_edge_has_zero_branch_length<P>(graph: &GraphAncestral, partitions: &[Arc<RwLock<P>>]) -> bool
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
///
/// Extracted from `run_optimize` so the Never-mode rejection path is unit
/// testable without standing up partitions and reading files.
pub(crate) fn apply_initial_guess_mode<P>(
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

fn write_graph<N, E, D>(outdir: impl AsRef<Path>, graph: &Graph<N, E, D>) -> Result<(), Report>
where
  N: GraphNode + NodeToNwk + Serialize,
  E: GraphEdge + EdgeToNwk + Serialize,
  D: Send + Sync + Default + Serialize,
{
  // json_write_file(
  //   outdir.as_ref().join("annotated_tree.graph.json"),
  //   &graph,
  //   JsonPretty(true),
  // )?;

  nwk_write_file(
    outdir.as_ref().join("annotated_tree.nwk"),
    graph,
    &NwkWriteOptions::default(),
  )?;

  nex_write_file(
    outdir.as_ref().join("annotated_tree.nexus"),
    graph,
    &NexWriteOptions::default(),
  )?;

  Ok(())
}
