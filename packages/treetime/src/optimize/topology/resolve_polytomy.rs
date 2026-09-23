use crate::make_internal_report;
use crate::optimize::params::TopologyOps;
use crate::optimize::topology::collapse::collapse_edge;
use crate::optimize::topology::hoist_reversions::{
  count_child_reversions, hoist_reverting_child, slide_bifurcating_root_for_child,
};
use crate::optimize::topology::merge_shared_mutations::merge_single_polytomy;
use crate::optimize::topology::polytomy_nodes::find_polytomy_nodes;
use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
use crate::partition::storage::sparse::SparseNodeState;
use eyre::Report;
use log::debug;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;

pub fn resolve_polytomies(
  graph: &mut Graph,
  sparse: &mut [PartitionMarginalSparse],
  node_states: &mut [BTreeMap<GraphNodeKey, SparseNodeState>],
  topology_ops: TopologyOps,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<usize, Report> {
  if sparse.is_empty() || !(topology_ops.merge_siblings || topology_ops.flip_parent_child) {
    return Ok(0);
  }

  let mut total_changed = 0;
  loop {
    let polytomy_keys = find_polytomy_nodes(graph);
    let mut round_changed = 0;
    for node_key in polytomy_keys {
      if resolve_one(graph, sparse, node_states, node_key, topology_ops, branch_lengths)? {
        round_changed += 1;
      }
    }
    if round_changed == 0 {
      break;
    }
    total_changed += round_changed;
  }

  if total_changed > 0 {
    debug!("Resolved {total_changed} polytomies via merge/hoist/retire");
  }

  Ok(total_changed)
}

fn resolve_one(
  graph: &mut Graph,
  sparse: &mut [PartitionMarginalSparse],
  node_states: &mut [BTreeMap<GraphNodeKey, SparseNodeState>],
  node_key: GraphNodeKey,
  topology_ops: TopologyOps,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<bool, Report> {
  let preexisting: BTreeSet<GraphNodeKey> = graph.get_nodes().map(|node| node.key()).collect();

  let mut any_changed = false;
  loop {
    let merged = topology_ops.merge_siblings && merge_single_polytomy(graph, sparse, node_key, branch_lengths)? > 0;
    let hoisted = topology_ops.flip_parent_child
      && try_hoist_reverting_child(graph, sparse, node_states, node_key, branch_lengths)?;
    let retired = retire_created_helpers(graph, sparse, &preexisting, branch_lengths)?;

    if !(merged || hoisted || retired) {
      break;
    }
    any_changed = true;
  }

  Ok(any_changed)
}

fn try_hoist_reverting_child(
  graph: &mut Graph,
  sparse: &mut [PartitionMarginalSparse],
  node_states: &mut [BTreeMap<GraphNodeKey, SparseNodeState>],
  node_key: GraphNodeKey,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<bool, Report> {
  let Some(parent_edge_key) = single_inbound_edge(graph, node_key) else {
    return Ok(false);
  };
  let degree_out = graph.get_node(node_key).map_or(0, |node| node.degree_out());
  if degree_out < 2 {
    return Ok(false);
  }
  let root_and_sibling = bifurcating_root_sibling_edge(graph, parent_edge_key)?;
  let sibling_edge_key = root_and_sibling.map(|(_, sibling_edge_key)| sibling_edge_key);
  let Some(child_edge_key) = best_reverting_child(graph, sparse, node_key, parent_edge_key, sibling_edge_key) else {
    return Ok(false);
  };
  if let Some((root_key, sibling_edge_key)) = root_and_sibling {
    slide_bifurcating_root_for_child(
      sparse,
      node_states,
      root_key,
      parent_edge_key,
      sibling_edge_key,
      child_edge_key,
    )?;
  }
  hoist_reverting_child(graph, sparse, parent_edge_key, child_edge_key, branch_lengths)?;
  Ok(true)
}

fn bifurcating_root_sibling_edge(
  graph: &Graph,
  parent_edge_key: GraphEdgeKey,
) -> Result<Option<(GraphNodeKey, GraphEdgeKey)>, Report> {
  let root_key = graph.get_source_node_key(parent_edge_key)?;
  let root = graph
    .get_node(root_key)
    .ok_or_else(|| make_internal_report!("Node {root_key} not found"))?;
  if !root.is_root() || root.degree_out() != 2 {
    return Ok(None);
  }
  let sibling_edge_key = root
    .outbound()
    .iter()
    .copied()
    .find(|&edge_key| edge_key != parent_edge_key);
  Ok(sibling_edge_key.map(|sibling_edge_key| (root_key, sibling_edge_key)))
}

fn best_reverting_child(
  graph: &Graph,
  sparse: &[PartitionMarginalSparse],
  node_key: GraphNodeKey,
  parent_edge_key: GraphEdgeKey,
  sibling_edge_key: Option<GraphEdgeKey>,
) -> Option<GraphEdgeKey> {
  let child_edges = graph.get_node(node_key)?.outbound().to_vec();

  let mut best: Option<(usize, GraphEdgeKey)> = None;
  for child_edge_key in child_edges {
    let reversions = count_child_reversions(sparse, parent_edge_key, sibling_edge_key, child_edge_key);
    if reversions == 0 {
      continue;
    }
    let better = match best {
      Some((best_reversions, best_key)) => {
        best_reversions > reversions || (best_reversions == reversions && best_key <= child_edge_key)
      },
      None => false,
    };
    if !better {
      best = Some((reversions, child_edge_key));
    }
  }

  best.map(|(_, key)| key)
}

fn retire_created_helpers(
  graph: &mut Graph,
  sparse: &mut [PartitionMarginalSparse],
  preexisting: &BTreeSet<GraphNodeKey>,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<bool, Report> {
  let mut retired = false;
  loop {
    let candidate = graph.get_edges().find_map(|edge| {
      let target_key = edge.target();
      if preexisting.contains(&target_key) {
        return None;
      }
      let target_is_leaf = graph.get_node(target_key).is_some_and(|node| node.is_leaf());
      if target_is_leaf {
        return None;
      }
      let edge_key = edge.key();
      let mutation_free = sparse.iter().all(|partition| match partition.obs_edges.get(&edge_key) {
        Some(edge_data) => edge_data.fitch_subs().is_empty() && edge_data.indels.is_empty(),
        None => true,
      });
      mutation_free.then_some(edge_key)
    });

    match candidate {
      Some(edge_key) => {
        collapse_edge(graph, sparse, edge_key, branch_lengths)?;
        retired = true;
      },
      None => break,
    }
  }
  Ok(retired)
}

fn single_inbound_edge(graph: &Graph, node_key: GraphNodeKey) -> Option<GraphEdgeKey> {
  let node = graph.get_node(node_key)?;
  match node.inbound() {
    [edge_key] => Some(*edge_key),
    _ => None,
  }
}
