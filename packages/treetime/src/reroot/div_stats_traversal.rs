use crate::make_report;
use crate::reroot::div_stats::DivStats;
use crate::reroot::traits::RootStats;
use crate::reroot::variance::VarianceModel;
use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::edge::{GraphEdge, GraphEdgeKey, HasBranchLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, GraphNodeKey};

/// Per-edge directional `DivStats` messages plus the aggregate at the current root.
pub struct DivStatsField {
  /// For each edge, `(to_parent, to_child)`: the child subtree as seen at the
  /// child node, and the rest of the tree as seen at the parent node.
  pub edge_stats: BTreeMap<GraphEdgeKey, (DivStats, DivStats)>,
  /// Aggregate statistics at the current root (the search baseline).
  pub root_stats: DivStats,
}

/// Compute `DivStats` messages for every edge by message passing over branch lengths.
///
/// Two passes mirror the clock regression: a leaves-to-root pass accumulates each
/// subtree message at its child node (`to_parent`) and propagates it across the
/// branch (`from_child`); a root-to-leaves pass derives the complementary
/// rest-of-tree message (`to_child`) by subtracting a child's contribution from
/// the node aggregate. Statistics are returned in maps rather than stored on the
/// graph, since the optimize payloads carry no message fields.
pub fn compute_div_stats<N, E, D>(graph: &Graph<N, E, D>, variance: &VarianceModel) -> Result<DivStatsField, Report>
where
  N: GraphNode,
  E: GraphEdge + HasBranchLength,
  D: Send + Sync,
{
  // Leaves-to-root order (children precede parents).
  let mut backward_order: Vec<GraphNodeKey> = Vec::new();
  graph.iter_breadth_first_backward(|n| {
    backward_order.push(n.key);
    Ok(())
  })?;

  let mut to_parent: BTreeMap<GraphEdgeKey, DivStats> = BTreeMap::new();
  let mut from_child: BTreeMap<GraphEdgeKey, DivStats> = BTreeMap::new();
  let mut root_stats = DivStats::default();

  // Backward pass.
  for &node_key in &backward_order {
    let NodeTopology { is_leaf, is_root, child_edges, parent_edge } = node_topology(graph, node_key)?;

    if is_root {
      root_stats = sum_children(&from_child, &child_edges)?;
      continue;
    }

    let parent_edge = parent_edge.ok_or_else(|| make_report!("Non-root node {node_key} has no parent edge"))?;
    let branch_length = branch_length_of(graph, parent_edge)?;

    if is_leaf {
      from_child.insert(
        parent_edge,
        DivStats::leaf(None, branch_length, variance.leaf_branch(branch_length)),
      );
      // Unused by the cost function (leaf edges take the leaf path), but recorded
      // so every edge has an entry: the leaf as seen at its own node, distance 0.
      to_parent.insert(parent_edge, DivStats::leaf(None, 0.0, variance.leaf_branch(0.0)));
    } else {
      let subtree = sum_children(&from_child, &child_edges)?;
      to_parent.insert(parent_edge, subtree);
      from_child.insert(
        parent_edge,
        subtree.propagate(branch_length, variance.branch(branch_length)),
      );
    }
  }

  // Forward pass (root-to-leaves): complementary rest-of-tree messages.
  let mut to_child: BTreeMap<GraphEdgeKey, DivStats> = BTreeMap::new();

  for &node_key in backward_order.iter().rev() {
    let NodeTopology { is_root, child_edges, parent_edge, .. } = node_topology(graph, node_key)?;

    // Full aggregate at this node: its own subtree plus the rest of the tree
    // propagated down the parent branch.
    let full = if is_root {
      root_stats
    } else {
      let parent_edge = parent_edge.ok_or_else(|| make_report!("Non-root node {node_key} has no parent edge"))?;
      let branch_length = branch_length_of(graph, parent_edge)?;
      let tp = *to_parent
        .get(&parent_edge)
        .ok_or_else(|| make_report!("Missing to_parent for edge {parent_edge}"))?;
      let tc = *to_child
        .get(&parent_edge)
        .ok_or_else(|| make_report!("Missing to_child for edge {parent_edge}"))?;
      tp + tc.propagate(branch_length, variance.branch(branch_length))
    };

    for child_edge in &child_edges {
      let fc = *from_child
        .get(child_edge)
        .ok_or_else(|| make_report!("Missing from_child for edge {child_edge}"))?;
      to_child.insert(*child_edge, full - fc);
    }
  }

  let edge_stats = to_parent
    .into_iter()
    .map(|(edge, tp)| {
      let tc = *to_child
        .get(&edge)
        .ok_or_else(|| make_report!("Missing to_child for edge {edge}"))?;
      Ok((edge, (tp, tc)))
    })
    .collect::<Result<BTreeMap<_, _>, Report>>()?;

  Ok(DivStatsField { edge_stats, root_stats })
}

fn sum_children(
  from_child: &BTreeMap<GraphEdgeKey, DivStats>,
  child_edges: &[GraphEdgeKey],
) -> Result<DivStats, Report> {
  let mut acc = DivStats::default();
  for edge in child_edges {
    acc = acc
      + *from_child
        .get(edge)
        .ok_or_else(|| make_report!("Missing child message for edge {edge}"))?;
  }
  Ok(acc)
}

struct NodeTopology {
  is_leaf: bool,
  is_root: bool,
  child_edges: Vec<GraphEdgeKey>,
  parent_edge: Option<GraphEdgeKey>,
}

fn node_topology<N, E, D>(graph: &Graph<N, E, D>, node_key: GraphNodeKey) -> Result<NodeTopology, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Send + Sync,
{
  let node = graph
    .get_node(node_key)
    .ok_or_else(|| make_report!("Node not found: {node_key}"))?;
  let node = node.read_arc();
  Ok(NodeTopology {
    is_leaf: node.is_leaf(),
    is_root: node.is_root(),
    child_edges: node.outbound().to_vec(),
    parent_edge: node.inbound().first().copied(),
  })
}

fn branch_length_of<N, E, D>(graph: &Graph<N, E, D>, edge_key: GraphEdgeKey) -> Result<f64, Report>
where
  N: GraphNode,
  E: GraphEdge + HasBranchLength,
  D: Send + Sync,
{
  graph
    .get_edge(edge_key)
    .ok_or_else(|| make_report!("Edge not found: {edge_key}"))?
    .read_arc()
    .payload()
    .read_arc()
    .branch_length()
    .ok_or_else(|| make_report!("Edge {edge_key} has no branch length"))
}
