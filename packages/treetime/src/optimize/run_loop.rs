use crate::ancestral::marginal::profile_branch_lengths;
use crate::ancestral::pipeline::{DenseReconstruction, SparseReconstruction};
use crate::gtr::gtr::GTR;
use crate::optimize::branch_length::invalid_branch_length_descriptions;
use crate::optimize::dispatch::initial_guess_mixed;
use crate::optimize::dispatch::run_optimize_mixed_inner;
use crate::optimize::gather::{gather_edge_contributions, gather_edge_indel_counts, total_sequence_length};
use crate::optimize::indel::{estimate_indel_rate, total_indel_log_lh};
use crate::optimize::iteration::apply_damping;
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

/// Run a marginal update over every sparse reconstruction, returning the updated reconstructions and
/// the summed substitution log likelihood.
pub fn marginal_update_sparse(
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

/// Run a marginal update over every dense reconstruction, returning the updated reconstructions and
/// the summed substitution log likelihood.
pub fn marginal_update_dense(
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
/// - Dense partitions have populated profiles ([`initialize_marginal`] + [`marginal_update`]).
/// - `mixed_partitions` contains both partition families (see [`collect_optimize_partitions`]).
/// - Each partition's `gtr` field is the final model (resolved before partition construction).
///
/// The indel rate is estimated once before the first iteration and held fixed
/// throughout the loop, removing the feedback path where branch-length changes
/// shift the rate estimate. When `no_indels` is true, the rate is zero and
/// the Poisson indel term drops out of both the likelihood evaluation and the
/// per-edge optimizer.
///
/// Branch lengths are the loop's source of truth, held in a map keyed by edge id. The initial
/// map is collected from the tree (the Newick-derived lengths after any initial guess and
/// reroot); the marginal reconstruction reads it directly rather than reaching back off the
/// graph edge, and the best map seen so far drives the worsening and numerical-failure rollback.
/// The graph edge branch length is written once, at loop end, so the output writers still reading
/// it see the final optimized (or rolled-back best) lengths.
///
/// Per-iteration sequence:
///
/// 1. Run `marginal_update` on sparse and dense partitions and sum the joint
///    substitution + indel log-likelihood using the pre-computed indel rate.
/// 2. Check three stopping conditions (converged, oscillating, worsened).
/// 3. Save the current branch-length map for damping.
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
  // Owned working copy of the routed names, refreshed from each `assign_node_names` the topology
  // cleanup runs and returned so `run_optimize` reads the final tree's names.
  let mut names = names.clone();
  // The loop's source of truth for branch lengths, keyed by edge id, holding `Option<f64>` so a
  // missing weight stays `None` end to end. Supplied by the caller (the parsed lengths after any
  // initial guess and reroot), then updated in place by the per-edge optimizer, damping, and topology
  // cleanup (topology producers insert new-edge keys and drop removed ones). The marginal
  // reconstruction reads the derived per-edge length (see [`profile_branch_lengths`]).
  let mut branch_lengths = branch_lengths;
  // The loop owns the partitions and returns them: every marginal update consumes the reconstructions
  // and hands back new ones, so no stage observes a half-updated reconstruction.
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
  // Best branch lengths seen so far, keyed by edge id. `None` until the first finite iteration
  // and reset whenever a topology change makes the keys stale, so a rollback never restores
  // lengths keyed to a superseded tree. On a worsening or numerically failed iteration the loop
  // recomputes the partitions from this map rather than restoring snapshotted partition state.
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
        let marginal_bl = profile_branch_lengths(&branch_lengths);
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
        let marginal_bl = profile_branch_lengths(&branch_lengths);
        (sparse_partitions, _) = marginal_update_sparse(graph, &marginal_bl, sparse_partitions)?;
        (dense_partitions, _) = marginal_update_dense(graph, &marginal_bl, dense_partitions)?;
      }
      stopped_at = Some((i, ConvergenceReason::Worsened));
      break;
    }

    // The per-edge optimizer, damping, and topology cleanup all read and update the
    // branch-length map in place; the map is the loop's source of truth.
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

  // The final branch-length map is the loop's result. On a normal exit it holds the optimized
  // lengths; after a rollback it holds the recovered best lengths. Callers read it directly (the
  // gather and the post-loop marginal pass).
  Ok(OptimizeLoopResult {
    sparse_partitions,
    dense_partitions,
    branch_lengths,
    names,
    lh_history,
    stopped_at,
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

/// Result of [`run_optimize_loop`].
#[derive(Clone, Debug, Default)]
pub struct OptimizeLoopResult {
  /// The sparse reconstructions the loop was given, at the state its final marginal update and
  /// topology cleanup left them in.
  pub sparse_partitions: Vec<SparseReconstruction>,

  /// The dense reconstructions the loop was given, at the state its final marginal update and
  /// topology cleanup left them in.
  pub dense_partitions: Vec<DenseReconstruction>,

  /// Final optimized (or rolled-back best) branch lengths, keyed by edge id. This is the loop's
  /// source of truth; `run_optimize` feeds it to the post-loop marginal pass and the output gather.
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,

  /// Final node names, keyed by node id, refreshed from each `assign_node_names` the loop's topology
  /// cleanup runs. `run_optimize` returns these.
  pub names: BTreeMap<GraphNodeKey, Option<String>>,

  /// Total log-likelihood recorded at the start of each iteration, before that iteration's
  /// branch-length update. Length equals the number of iterations actually executed
  /// (including the final iteration that triggered the convergence break, if any).
  #[allow(dead_code, reason = "read only in test code")]
  pub lh_history: Vec<LogLh>,

  /// Iteration index (0-based) and reason the loop stopped early.
  /// `None` if the loop exhausted `max_iter` without meeting any stopping criterion.
  #[allow(dead_code, reason = "read only in test code")]
  pub stopped_at: Option<(usize, ConvergenceReason)>,
}

#[derive(Clone, Copy, Debug, Default)]
struct OptimizeIterationLikelihood {
  sparse_lh: LogLh,
  dense_lh: LogLh,
  indel_lh: LogLh,
  total_lh: LogLh,
  indel_rate: f64,
}

/// Run one iteration's marginal update over both partition families and score the result.
fn compute_iteration(
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  sparse_partitions: Vec<SparseReconstruction>,
  dense_partitions: Vec<DenseReconstruction>,
  indel_rate: f64,
  no_indels: bool,
) -> Result<OptimizeIteration, Report> {
  let marginal_bl = profile_branch_lengths(branch_lengths);
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

/// One iteration's updated reconstructions together with the likelihood terms they produced.
struct OptimizeIteration {
  sparse_partitions: Vec<SparseReconstruction>,
  dense_partitions: Vec<DenseReconstruction>,
  likelihood: OptimizeIterationLikelihood,
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
  graph: &Graph,
  sparse_partitions: &[SparseReconstruction],
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Vec<GraphEdgeKey> {
  graph
    .get_edges()
    .iter()
    .filter_map(|edge_ref| {
      let edge = edge_ref.read_arc();
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

/// Collapse zero-optimal internal edges, then resolve polytomies in the resulting tree.
///
/// Two topology-cleanup steps run per iteration:
///
/// 1. Collapse zero-optimal internal edges. Analogous to v0's `prune_short_branches()`:
///    edges whose optimal branch length is zero are collapsed, simplifying the tree
///    progressively across iterations and producing polytomies.
/// 2. Resolve polytomies ([`resolve_polytomies`], sparse partitions only): merge siblings
///    that share substitutions, hoist a reverting child under a new node to remove one
///    mutation per reversion, and retire the helper nodes left behind.
///
/// Step 2 runs every iteration, not only after a collapse fired, because reverting-child and
/// shared-mutation polytomies exist in the input and in polytomies formed by earlier
/// iterations, independent of any zero-optimal collapse in the current one. The combined
/// (mutation count, node count) potential is monotone non-increasing, so the two steps cannot
/// oscillate across iterations even though the hoist deliberately relocates mutation-carrying
/// edges that [`find_zero_optimal_internal_edges`] refuses to collapse.
///
/// Returns true if any topology change occurred.
pub fn prune_and_merge_in_loop(
  graph: &mut Graph,
  sparse_partitions: Vec<SparseReconstruction>,
  dense_partitions: Vec<DenseReconstruction>,
  zero_optimal_edges: &[GraphEdgeKey],
  topology_ops: TopologyOps,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  names: &mut BTreeMap<GraphNodeKey, Option<String>>,
) -> Result<TopologyCleanup, Report> {
  // The topology moves read and rewrite the sparse observations, and one of them (the
  // bifurcating-root slide) keeps the root node state in step with them, so each reconstruction is
  // split here into the values a structural change carries across -- observations and node states --
  // and the per-edge results of the last update, which it does not. With no change both parts go back
  // together unchanged. With a change the node states are reconciled to the new node set and the
  // per-edge results are left behind, because they describe the superseded topology.
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
    // Override damped branch lengths back to zero for edges the optimizer identified
    // as zero-optimal. Damping is a convergence aid for continuous values; it should
    // not prevent collapsing degenerate edges.
    for &edge_key in zero_optimal_edges {
      if branch_lengths.contains_key(&edge_key) {
        branch_lengths.insert(edge_key, Some(0.0));
      }
    }

    let mut collapsed = 0_usize;
    for &edge_key in zero_optimal_edges {
      // Edge may already be gone if a prior collapse in this batch removed it
      // (e.g., collapsing a parent also removed a child edge that was zero-optimal)
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
  // One complete remap at the structural-operation boundary: the topology moves kept the sparse
  // observations precise, but the node-state maps still key to the pre-change topology. Reconcile them
  // to the current node set (dropping removed nodes, seeding placeholders for created ones); the next
  // marginal update rebuilds a complete set of per-edge results over the new topology.
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

/// The reconstructions a topology cleanup returns, together with whether it changed the topology.
pub struct TopologyCleanup {
  pub sparse_partitions: Vec<SparseReconstruction>,
  pub dense_partitions: Vec<DenseReconstruction>,
  pub topology_changed: bool,
}

/// Put the split sparse reconstructions back together unchanged, for the case where the topology
/// moves changed nothing.
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

/// Put the split dense reconstructions back together unchanged, for the case where the topology moves
/// changed nothing.
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
pub fn any_indel_edge_has_zero_branch_length(
  graph: &Graph,
  indel_counts: &BTreeMap<GraphEdgeKey, usize>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
) -> bool {
  graph.get_edges().iter().any(|edge_ref| {
    let edge = edge_ref.read_arc();
    let edge_key = edge.key();
    let bl = branch_lengths[&edge_key].unwrap_or(0.0);
    if bl != 0.0 {
      return false;
    }
    indel_counts[&edge_key] > 0
  })
}

/// Apply the initial-guess mode to the graph.
///
/// - `Auto`: estimate only edges with missing or invalid branch lengths.
/// - `Always`: estimate all edges, overwriting existing values.
/// - `Never`: keep input branch lengths; error if any edge lacks a usable
///   value (None/NaN) or has zero branch length while carrying indels (the
///   Poisson indel log-likelihood diverges at $t = 0$ when $k > 0$).
#[allow(clippy::too_many_arguments)]
pub fn apply_initial_guess_mode(
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
      false,
      no_indels,
      branch_lengths,
    ),
    InitialGuessMode::Always => initial_guess_mixed(
      graph,
      total_length,
      indel_counts,
      sub_counts,
      effective_lengths,
      true,
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
      // `Never` makes no branch-length writes, so the input tree's lengths are the ones checked.
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
pub fn normalize_partition_rates(
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
