use crate::optimize::args::InitialGuessMode;
use crate::optimize::optimize_unified::initial_guess_mixed;
use crate::representation::algo::topology_cleanup::collapse::collapse_edge;
use crate::representation::algo::topology_cleanup::merge_shared_mutations::merge_shared_mutation_branches;
use crate::representation::partition::marginal_dense::PartitionMarginalDense;
use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
use crate::representation::partition::traits::{HasGtr, PartitionOptimizeOps};
use crate::representation::payload::ancestral::GraphAncestral;
use eyre::Report;
use itertools::{Itertools, izip};
use log::debug;
use num_traits::pow::pow;
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_graph::edge::{GraphEdge, GraphEdgeKey, HasBranchLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNode;
use treetime_utils::make_error;

/// Save current branch lengths for all edges in graph traversal order.
pub fn save_branch_lengths<N, E, D>(graph: &Graph<N, E, D>) -> Vec<f64>
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
/// from the same graph (same edge count and order).
pub fn restore_branch_lengths<N, E, D>(graph: &Graph<N, E, D>, saved: &[f64])
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
pub const DAMPING_FLOOR: f64 = 0.01;

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
pub fn apply_damping<N, E, D>(graph: &Graph<N, E, D>, old_branch_lengths: &[f64], damping: f64, iteration: usize)
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
