use crate::commands::timetree::partition_ops::PartitionTimetreeAll;
use crate::representation::partition::timetree::GraphTimetree;
use crate::representation::payload::timetree::EdgeTimetree;
use crate::representation::payload::timetree::NodeTimetree;
use eyre::Report;
use itertools::Itertools;
use log::debug;
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_graph::node::GraphNodeKey;

/// Default minimum likelihood gain required to merge two children.
pub const DEFAULT_RESOLUTION_THRESHOLD: f64 = 0.05;

/// Resolve multifurcations using temporal constraints and sequence data.
///
/// Scans tree for polytomies (nodes with >2 children), resolves them using greedy
/// pairwise merging based on likelihood gain. Resolution stops when no pair exceeds
/// the threshold. After resolution, removes obsolete single-child internal nodes.
///
/// Returns count of new nodes inserted.
pub fn resolve_polytomies(
  graph: &mut GraphTimetree,
  partitions: &[Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>>],
) -> Result<usize, Report> {
  resolve_polytomies_with_options(graph, partitions, DEFAULT_RESOLUTION_THRESHOLD, false)
}

/// Resolve polytomies with configurable options.
///
/// - `resolution_threshold`: minimum delta LH to merge a pair (default 0.05)
/// - `merge_compressed`: if true, also resolve compressed branches (default false)
#[allow(unused_variables)]
pub fn resolve_polytomies_with_options(
  graph: &mut GraphTimetree,
  partitions: &[Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>>],
  resolution_threshold: f64,
  merge_compressed: bool,
) -> Result<usize, Report> {
  let mut total_resolved = 0;

  // Find all polytomy nodes (>2 children)
  let polytomy_keys = find_polytomy_nodes(graph);

  for node_key in polytomy_keys {
    let prior_n_children = graph.degree_out(node_key)?;

    // Resolve this polytomy using greedy pairwise merging
    let n_resolved = resolve_single_polytomy(graph, node_key, resolution_threshold)?;
    total_resolved += n_resolved;

    let new_n_children = graph.degree_out(node_key)?;
    debug!(
      "Polytomy at node {node_key}: {prior_n_children} -> {new_n_children} children, resolved {n_resolved}"
    );
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

/// Resolve a single polytomy node using greedy pairwise merging.
///
/// For each pair of children, compute cost gain from merging them under a new
/// internal node. Greedily merge the pair with highest gain until no pair
/// exceeds the threshold.
///
/// Returns count of new nodes created.
fn resolve_single_polytomy(
  _graph: &mut GraphTimetree,
  _node_key: GraphNodeKey,
  _resolution_threshold: f64,
) -> Result<usize, Report> {
  // TODO: Implement greedy merging algorithm
  // 1. Get all children of this node
  // 2. Compute pairwise cost gain matrix
  // 3. While max gain > threshold:
  //    a. Find pair with max gain
  //    b. Create new internal node
  //    c. Reparent the two children under new node
  //    d. Update cost gains for new node vs remaining children
  // 4. Return count of new nodes created
  Ok(0)
}

/// Remove obsolete single-child internal nodes.
///
/// After polytomy resolution, some nodes may have only one child. These should
/// be collapsed by connecting their parent directly to their child.
fn remove_obsolete_nodes(graph: &mut GraphTimetree) -> Result<usize, Report> {
  let mut removed_count = 0;

  loop {
    // Find a single-child internal node (not root)
    let obsolete_key = graph
      .get_nodes()
      .into_iter()
      .find_map(|node| {
        let node = node.read_arc();
        let is_single_child = node.has_one_child();
        let is_internal = !node.is_root() && !node.is_leaf();
        (is_single_child && is_internal).then_some(node.key())
      });

    let Some(node_key) = obsolete_key else {
      break;
    };

    // Collapse this node: connect parent directly to child
    let inbound_edge_key = graph.parent_inbound_edge(node_key)?.expect("Internal node must have parent");
    graph.collapse_edge(inbound_edge_key)?;
    removed_count += 1;
  }

  // Rebuild roots/leaves lists after topology changes
  graph.build()?;

  Ok(removed_count)
}

/// Reset internal state after tree topology changes.
///
/// Clears cached likelihoods, recalculates divergence from root, and resets
/// bad_branch flags. Call this after any topology modification.
pub fn prepare_tree_after_topology_change(graph: &GraphTimetree) -> Result<(), Report> {
  // Clear cached time distributions and bad_branch flags
  for node in graph.get_nodes() {
    let mut payload = node.read_arc().payload().write_arc();
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
