use crate::representation::partition::timetree::GraphTimetree;
use crate::representation::partition::traits::PartitionTimetreeAll;
use crate::representation::payload::timetree::EdgeTimetree;
use crate::representation::payload::timetree::NodeTimetree;
use argmin::core::{CostFunction, Executor};
use argmin::solver::brent::BrentOpt;
use eyre::Report;
use itertools::Itertools;
use log::debug;
use ndarray::Array2;
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_distribution::Distribution;
use treetime_graph::edge::{GraphEdgeKey, HasBranchLength};
use treetime_graph::node::GraphNodeKey;

/// Default minimum likelihood gain required to merge two children.
pub const DEFAULT_RESOLUTION_THRESHOLD: f64 = 0.05;

/// Resolve multifurcations using temporal constraints and sequence data.
///
/// Scans tree for polytomies (nodes with >2 children), resolves them using greedy
/// pairwise merging based on likelihood gain. Only "stretched" children (fewer
/// mutations than clock predicts) are merged by default, matching v0 behavior.
/// Resolution stops when no pair exceeds the threshold. After resolution, removes
/// obsolete single-child internal nodes.
///
/// `zero_branch_slope` = `clock_rate * sequence_length`.
/// `clock_rate` = substitutions per site per time unit.
pub fn resolve_polytomies(
  graph: &mut GraphTimetree,
  partitions: &[Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>>],
  zero_branch_slope: f64,
  clock_rate: f64,
) -> Result<usize, Report> {
  resolve_polytomies_with_options(graph, partitions, DEFAULT_RESOLUTION_THRESHOLD, zero_branch_slope, clock_rate, false)
}

/// Resolve polytomies with configurable options.
///
/// - `resolution_threshold`: minimum delta LH to merge a pair (default 0.05)
/// - `zero_branch_slope`: expected substitutions per unit time (`clock_rate * sequence_length`)
/// - `clock_rate`: substitutions per site per time unit, used to classify children
///   as "stretched" (mutation_length < clock_length) vs "compressed"
/// - `merge_compressed`: if true, also merge compressed children after stretched
pub fn resolve_polytomies_with_options(
  graph: &mut GraphTimetree,
  partitions: &[Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>>],
  resolution_threshold: f64,
  zero_branch_slope: f64,
  clock_rate: f64,
  merge_compressed: bool,
) -> Result<usize, Report> {
  let mut total_resolved = 0;

  // Find all polytomy nodes (>2 children)
  let polytomy_keys = find_polytomy_nodes(graph);

  for node_key in polytomy_keys {
    let prior_n_children = graph.degree_out(node_key)?;

    let n_resolved = resolve_single_polytomy(
      graph,
      node_key,
      resolution_threshold,
      zero_branch_slope,
      clock_rate,
      merge_compressed,
    )?;
    total_resolved += n_resolved;

    let new_n_children = graph.degree_out(node_key)?;
    debug!("Polytomy at node {node_key}: {prior_n_children} -> {new_n_children} children, resolved {n_resolved}");
  }

  // Remove obsolete nodes (single-child internal nodes created during resolution)
  let obsolete_count = remove_obsolete_nodes(graph)?;
  if obsolete_count > 0 {
    debug!("Removed {obsolete_count} obsolete single-child nodes");
  }

  if total_resolved > 0 {
    debug!("Polytomy resolution introduced {total_resolved} new nodes");
  } else {
    debug!("No polytomies to resolve");
  }

  Ok(total_resolved)
}

/// Find all nodes with more than 2 children (polytomies).
fn find_polytomy_nodes(graph: &GraphTimetree) -> Vec<GraphNodeKey> {
  graph
    .get_nodes()
    .into_iter()
    .filter_map(|node| {
      let node = node.read_arc();
      (node.degree_out() > 2).then_some(node.key())
    })
    .collect_vec()
}

/// Child info needed for polytomy resolution.
struct ChildInfo {
  node_key: GraphNodeKey,
  edge_key: GraphEdgeKey,
  time: f64,
  branch_dist: Arc<Distribution>,
  /// Observed substitutions per site on this branch.
  mutation_length: f64,
  /// Clock-expected substitutions per site (`time_length * clock_rate`).
  clock_length: f64,
}

impl ChildInfo {
  /// A "stretched" branch has fewer mutations than the clock predicts.
  fn is_stretched(&self) -> bool {
    self.mutation_length < self.clock_length
  }
}

/// Result of computing cost gain for merging two children.
struct MergeCandidate {
  optimal_time: f64,
  cost_gain: f64,
}

/// Result of finding the best pair of children to merge.
struct BestPair {
  child_idx_a: usize,
  child_idx_b: usize,
  likelihood_gain: f64,
}

/// Resolve a single polytomy node using greedy pairwise merging.
///
/// Children are classified as "stretched" (fewer mutations than clock predicts)
/// or "compressed". Only stretched children are merged by default. When
/// `merge_compressed` is true, compressed children are also merged in a second pass.
///
/// Within each group, the loop merges the pair with the highest likelihood gain
/// until gain drops below the threshold or fewer than 2 candidates remain. When
/// merging a subset (not all children), the group is merged down to 1 remaining
/// child (creating one subtree node). When merging all children, resolution stops
/// at 2 remaining (normal bifurcation).
fn resolve_single_polytomy(
  graph: &mut GraphTimetree,
  node_key: GraphNodeKey,
  resolution_threshold: f64,
  zero_branch_slope: f64,
  clock_rate: f64,
  merge_compressed: bool,
) -> Result<usize, Report> {
  let mut nodes_created = 0;

  // Merge stretched children
  nodes_created += merge_child_group(graph, node_key, resolution_threshold, zero_branch_slope, clock_rate, true)?;

  // Optionally merge compressed children
  if merge_compressed {
    nodes_created +=
      merge_child_group(graph, node_key, resolution_threshold, zero_branch_slope, clock_rate, false)?;
  }

  Ok(nodes_created)
}

/// Merge one group (stretched or compressed) of a polytomy's children.
///
/// When `select_stretched` is true, only stretched children are candidates.
/// When false, only compressed children are candidates.
///
/// The stopping condition depends on whether the group covers all children:
/// - All children in group: stop at 2 remaining (normal bifurcation)
/// - Subset: stop at 1 remaining (merge entire subset into one subtree)
fn merge_child_group(
  graph: &mut GraphTimetree,
  node_key: GraphNodeKey,
  resolution_threshold: f64,
  zero_branch_slope: f64,
  clock_rate: f64,
  select_stretched: bool,
) -> Result<usize, Report> {
  let mut nodes_created = 0;

  loop {
    let parent_time = {
      let node = graph.get_node(node_key).expect("Node must exist");
      node.read_arc().payload().read_arc().time.unwrap_or(0.0)
    };

    let all_children = collect_children_info(graph, node_key, clock_rate)?;
    if all_children.len() <= 2 {
      break;
    }

    let candidates: Vec<_> = all_children
      .into_iter()
      .filter(|c| c.is_stretched() == select_stretched)
      .collect();

    if candidates.len() < 2 {
      break;
    }

    let n = candidates.len();
    let mut gains = Array2::<f64>::from_elem((n, n), f64::NEG_INFINITY);
    let mut optimal_times = Array2::<f64>::zeros((n, n));

    for i in 0..n {
      for j in (i + 1)..n {
        if let Some(candidate) = compute_merge_gain(&candidates[i], &candidates[j], parent_time, zero_branch_slope) {
          gains[[i, j]] = candidate.cost_gain;
          gains[[j, i]] = candidate.cost_gain;
          optimal_times[[i, j]] = candidate.optimal_time;
          optimal_times[[j, i]] = candidate.optimal_time;
        }
      }
    }

    let BestPair {
      child_idx_a,
      child_idx_b,
      likelihood_gain,
    } = find_best_pair(&gains);
    if likelihood_gain < resolution_threshold {
      debug!("Best gain {likelihood_gain:.4} below threshold {resolution_threshold:.4}, stopping");
      break;
    }

    let optimal_time = optimal_times[[child_idx_a, child_idx_b]];
    debug!(
      "Merging children {} and {} with gain {likelihood_gain:.4} at time {optimal_time:.4}",
      candidates[child_idx_a].node_key, candidates[child_idx_b].node_key
    );

    merge_children(
      graph,
      node_key,
      &candidates[child_idx_a],
      &candidates[child_idx_b],
      optimal_time,
      parent_time,
    )?;
    nodes_created += 1;
  }

  Ok(nodes_created)
}

/// Collect information about all children of a node.
fn collect_children_info(
  graph: &GraphTimetree,
  node_key: GraphNodeKey,
  clock_rate: f64,
) -> Result<Vec<ChildInfo>, Report> {
  let node = graph.get_node(node_key).expect("Node must exist");
  let node = node.read_arc();

  let mut children = Vec::new();
  for &edge_key in node.outbound() {
    let edge = graph.get_edge(edge_key).expect("Edge must exist");
    let edge = edge.read_arc();
    let child_key = edge.target();

    let child_node = graph.get_node(child_key).expect("Child must exist");
    let child_payload = child_node.read_arc().payload().read_arc();
    let edge_payload = edge.payload().read_arc();

    let time = child_payload.time.unwrap_or(0.0);
    let branch_dist = edge_payload
      .branch_length_distribution
      .clone()
      .unwrap_or_else(|| Arc::new(Distribution::Empty));
    let mutation_length = edge_payload.branch_length().unwrap_or(0.0);
    let time_length = edge_payload.time_length.unwrap_or(0.0);
    let clock_length = time_length * clock_rate;

    children.push(ChildInfo {
      node_key: child_key,
      edge_key,
      time,
      branch_dist,
      mutation_length,
      clock_length,
    });
  }

  Ok(children)
}

/// Compute cost gain from merging two children under a new internal node.
///
/// Cost gain = log(P_new / P_old) where:
/// - P_old = P(branch1_old) * P(branch2_old)
/// - P_new = P(branch1_new) * P(branch2_new) * P(new_branch)
///
/// The new branch has zero mutations, so its probability depends on time only.
///
/// Time convention: v1 uses calendar time where parent_time < child_time
/// (parent is in the past with smaller calendar year, children are in the future).
fn compute_merge_gain(
  child1: &ChildInfo,
  child2: &ChildInfo,
  parent_time: f64,
  zero_branch_slope: f64,
) -> Option<MergeCandidate> {
  // In calendar time: parent is older (smaller), children are younger (larger)
  // New node must be between parent and children: parent_time < new_time < min(child_times)
  let child_min_time = f64::min(child1.time, child2.time);
  if parent_time >= child_min_time {
    return None; // Invalid time constraints
  }

  // Define cost function for optimization
  // We minimize negative gain (i.e., maximize gain)
  let cost_fn = MergeCostFunction {
    child1_dist: &child1.branch_dist,
    child2_dist: &child2.branch_dist,
    child1_time: child1.time,
    child2_time: child2.time,
    parent_time,
    zero_branch_slope,
  };

  // Use Brent optimization to find optimal merge time
  // Bounds: parent_time < new_node_time < child_min_time
  let solver = BrentOpt::new(parent_time + 1e-10, child_min_time - 1e-10);
  let result = Executor::new(cost_fn.clone(), solver)
    .configure(|cfg| cfg.max_iters(50))
    .run();

  if let Ok(res) = result {
    let optimal_time = res
      .state()
      .best_param
      .unwrap_or_else(|| f64::midpoint(parent_time, child_min_time));
    let cost_gain = -res.state().best_cost; // We minimized negative gain
    Some(MergeCandidate {
      optimal_time,
      cost_gain,
    })
  } else {
    // Fallback: evaluate at midpoint
    let mid_time = f64::midpoint(parent_time, child_min_time);
    let cost_gain = -cost_fn.cost(&mid_time).unwrap_or(0.0);
    Some(MergeCandidate {
      optimal_time: mid_time,
      cost_gain,
    })
  }
}

/// Cost function for Brent optimization of merge time.
///
/// `zero_branch_slope` = `clock_rate * sequence_length`: expected substitutions per
/// unit time for the full alignment, scaling the penalty for zero-mutation branches.
#[derive(Clone)]
struct MergeCostFunction<'a> {
  child1_dist: &'a Distribution,
  child2_dist: &'a Distribution,
  child1_time: f64,
  child2_time: f64,
  parent_time: f64,
  zero_branch_slope: f64,
}

impl CostFunction for MergeCostFunction<'_> {
  type Param = f64;
  type Output = f64;

  fn cost(&self, t: &Self::Param) -> Result<Self::Output, argmin::core::Error> {
    let new_node_time = *t;

    // Branch lengths in calendar time: child_time - parent_time (positive since child > parent)
    // Old branches: from original parent directly to children
    let old_branch1 = self.child1_time - self.parent_time;
    let old_branch2 = self.child2_time - self.parent_time;

    // New branches: from new internal node to children
    let new_branch1 = self.child1_time - new_node_time;
    let new_branch2 = self.child2_time - new_node_time;
    // Branch from parent to new node
    let new_branch_to_parent = new_node_time - self.parent_time;

    // Evaluate log probabilities using branch length distributions
    // Distribution::eval returns probability, so we take log
    let log_prob_old1 = self.child1_dist.eval(old_branch1).unwrap_or(1e-10).ln();
    let log_prob_old2 = self.child2_dist.eval(old_branch2).unwrap_or(1e-10).ln();
    let log_prob_new1 = self.child1_dist.eval(new_branch1).unwrap_or(1e-10).ln();
    let log_prob_new2 = self.child2_dist.eval(new_branch2).unwrap_or(1e-10).ln();

    // The new branch from parent to new_node has zero mutations.
    // Under a Poisson model: P(0 mutations | mu, L, dt) = exp(-mu*L*dt),
    // so -log P(0) = mu*L*dt = zero_branch_slope * dt.
    let zero_branch_penalty = self.zero_branch_slope * new_branch_to_parent;

    // Cost gain = (new log prob) - (old log prob)
    // = (log_prob_new1 + log_prob_new2 - zero_branch_penalty) - (log_prob_old1 + log_prob_old2)
    let cost_gain = (log_prob_new1 + log_prob_new2 - zero_branch_penalty) - (log_prob_old1 + log_prob_old2);

    // Return negative gain for minimization
    Ok(-cost_gain)
  }
}

/// Find the pair with maximum cost gain.
fn find_best_pair(gains: &Array2<f64>) -> BestPair {
  let n = gains.nrows();
  let mut child_idx_a = 0;
  let mut child_idx_b = 1;
  let mut likelihood_gain = f64::NEG_INFINITY;

  for i in 0..n {
    for j in (i + 1)..n {
      if gains[[i, j]] > likelihood_gain {
        likelihood_gain = gains[[i, j]];
        child_idx_a = i;
        child_idx_b = j;
      }
    }
  }

  BestPair {
    child_idx_a,
    child_idx_b,
    likelihood_gain,
  }
}

/// Merge two children under a new internal node.
fn merge_children(
  graph: &mut GraphTimetree,
  parent_key: GraphNodeKey,
  child1: &ChildInfo,
  child2: &ChildInfo,
  new_node_time: f64,
  parent_time: f64,
) -> Result<(), Report> {
  // Create new internal node
  let new_payload = NodeTimetree {
    time: Some(new_node_time),
    ..NodeTimetree::default()
  };
  let new_node_key = graph.add_node(new_payload);

  // Create edge from parent to new node
  // Branch length in calendar time: child - parent (new_node is the "child" of parent)
  let new_branch_length = new_node_time - parent_time;
  let new_edge_payload = EdgeTimetree {
    time_length: Some(new_branch_length),
    ..EdgeTimetree::default()
  };
  graph.add_edge(parent_key, new_node_key, new_edge_payload)?;

  // Reparent child1: remove old edge from parent, create edge from new node
  // Branch length: child1.time - new_node_time
  let new_branch1_length = child1.time - new_node_time;
  graph.remove_edge(child1.edge_key)?;
  let child1_edge_payload = EdgeTimetree {
    time_length: Some(new_branch1_length),
    ..EdgeTimetree::default()
  };
  graph.add_edge(new_node_key, child1.node_key, child1_edge_payload)?;

  // Reparent child2: remove old edge from parent, create edge from new node
  // Branch length: child2.time - new_node_time
  let new_branch2_length = child2.time - new_node_time;
  graph.remove_edge(child2.edge_key)?;
  let child2_edge_payload = EdgeTimetree {
    time_length: Some(new_branch2_length),
    ..EdgeTimetree::default()
  };
  graph.add_edge(new_node_key, child2.node_key, child2_edge_payload)?;

  Ok(())
}

/// Remove obsolete single-child internal nodes.
///
/// After polytomy resolution, some nodes may have only one child. These should
/// be collapsed by connecting their parent directly to their child.
fn remove_obsolete_nodes(graph: &mut GraphTimetree) -> Result<usize, Report> {
  let mut removed_count = 0;

  loop {
    // Find a single-child internal node (not root)
    let obsolete_key = graph.get_nodes().into_iter().find_map(|node| {
      let node = node.read_arc();
      let is_single_child = node.has_one_child();
      let is_internal = !node.is_root() && !node.is_leaf();
      (is_single_child && is_internal).then_some(node.key())
    });

    let Some(node_key) = obsolete_key else {
      break;
    };

    // Collapse this node: connect parent directly to child
    let inbound_edge_key = graph
      .parent_inbound_edge(node_key)?
      .expect("Internal node must have parent");
    graph.collapse_edge(inbound_edge_key)?;
    removed_count += 1;
  }

  // Rebuild roots/leaves lists after topology changes
  graph.build()?;

  Ok(removed_count)
}

/// Reset internal state after tree topology changes.
///
/// Clears cached time distributions and bad_branch flags on internal nodes only.
/// Leaf nodes retain their date constraints (`time_distribution` from `load_date_constraints()`)
/// and exclusion flags (`bad_branch` from outlier detection or dateless leaves).
/// Edge distributions and messages are cleared unconditionally.
pub fn prepare_tree_after_topology_change(graph: &GraphTimetree) -> Result<(), Report> {
  // Clear cached time distributions and bad_branch flags on internal nodes only.
  // Leaf nodes must retain their date constraints and exclusion flags because
  // load_date_constraints() runs once at initialization and is not re-invoked.
  for node in graph.get_nodes() {
    let node = node.read_arc();
    if node.is_leaf() {
      continue;
    }
    let mut payload = node.payload().write_arc();
    payload.time_distribution = None;
    payload.bad_branch = false;
  }

  // Clear cached branch length distributions
  for edge in graph.get_edges() {
    let mut payload = edge.read_arc().payload().write_arc();
    payload.branch_length_distribution = None;
    payload.msg_to_parent = None;
  }

  Ok(())
}
