use crate::gtr::gtr::GTR;
use crate::partition::marginal::shared::data::IndexedMarginalPartition;
use crate::partition::marginal::shared::normalize::{
  forward_log_lh_add_normalization, forward_log_lh_remove_child, normalize_from_log, normalize_inplace,
};
use crate::partition::storage::dense::{DenseEdgePartition, DenseNodePartition, DenseSeqDistribution};
use eyre::Report;
use itertools::Itertools;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::pass::{GraphPass, GraphPassBackwardContext, GraphPassForwardContext, GraphPassNodeOutput};
use treetime_primitives::LogLh;

pub fn marginal_process_backward_indexed(
  partition: &mut dyn IndexedMarginalPartition,
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
) -> Result<(), Report> {
  let mut missing_nodes = BTreeMap::new();
  for node in graph.get_nodes() {
    let key = node.read_arc().key();
    if !partition.marginal_data().nodes.contains_key(&key) {
      missing_nodes.insert(key, partition.indexed_missing_node(key)?);
    }
  }
  partition.marginal_data_mut().nodes.append(&mut missing_nodes);
  let gtr = partition.marginal_data().gtr.clone();
  let min_branch_length = partition.marginal_data().min_branch_length;
  let pass = GraphPass::new(graph)?;
  let nodes = &partition.marginal_data().nodes;
  let edges = &partition.marginal_data().edges;
  let outputs = pass.map_backward(
    nodes,
    edges,
    |_| unreachable!("Missing nodes were initialized"),
    |context| marginal_process_node_backward_indexed(partition, &gtr, min_branch_length, branch_lengths, &context),
  )?;
  partition.marginal_data_mut().nodes = outputs.nodes;
  partition.marginal_data_mut().edges = outputs.edges;
  Ok(())
}

fn marginal_process_node_backward_indexed(
  partition: &dyn IndexedMarginalPartition,
  gtr: &GTR,
  min_branch_length: f64,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  context: &GraphPassBackwardContext<
    '_,
    DenseNodePartition,
    DenseEdgePartition,
    DenseNodePartition,
    DenseEdgePartition,
  >,
) -> Result<GraphPassNodeOutput<DenseNodePartition, DenseEdgePartition>, Report> {
  let mut node = context.input.clone();
  let msg_to_parent = if context.is_leaf {
    partition.indexed_leaf_profile(&node)?
  } else {
    // Children arrive in the graph's canonical `children_of` order, so the per-child log-space product
    // folds in the same order regardless of thread count, keeping the result byte-for-byte identical.
    let children = context.children.iter().map(|child| child.node).collect_vec();
    node.seq = partition.indexed_backward_internal(&children)?;

    let child_edges = context
      .children
      .iter()
      .map(|child| {
        child
          .edge
          .expect("Backward child edge message must be published before its parent")
      })
      .collect_vec();
    let first_edge = child_edges.first().expect("Internal node must have children");
    let mut log_dis = first_edge.msg_from_child.dis.mapv(f64::ln);
    for edge in child_edges.iter().skip(1) {
      log_dis += &edge.msg_from_child.dis.mapv(f64::ln);
    }
    let (dis, delta_ll) = normalize_from_log(&log_dis);
    let log_lh = child_edges.iter().map(|edge| edge.msg_from_child.log_lh).sum::<LogLh>() + LogLh::new(delta_ll);
    node.profile = DenseSeqDistribution {
      dis: dis.clone(),
      log_lh,
    };
    DenseSeqDistribution { dis, log_lh }
  };

  let parent_message = if context.is_root {
    let mut dis = &msg_to_parent.dis * &gtr.pi;
    let delta_ll = normalize_inplace(&mut dis);
    node.profile = DenseSeqDistribution {
      dis,
      log_lh: msg_to_parent.log_lh + LogLh::new(delta_ll),
    };
    None
  } else {
    // Clone this node's parent-edge input and overwrite only the two message fields, exactly as before.
    // The edge's other fields (indels, transmission, msg_to_child) carry forward-pass state across
    // timetree iterations and must survive the backward pass unchanged.
    let (edge_key, edge) = context.parent_edge.expect("Non-root node must own its parent edge");
    let mut edge = edge.clone();
    let branch_length = branch_lengths[&edge_key].max(min_branch_length);
    edge.msg_from_child = DenseSeqDistribution {
      dis: gtr.propagate_profile(&msg_to_parent.dis, branch_length, false),
      log_lh: msg_to_parent.log_lh,
    };
    edge.msg_to_parent = msg_to_parent;
    Some(edge)
  };

  Ok(GraphPassNodeOutput { node, parent_message })
}

pub fn marginal_process_forward_indexed(
  partition: &mut dyn IndexedMarginalPartition,
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
) -> Result<(), Report> {
  let gtr = partition.marginal_data().gtr.clone();
  let min_branch_length = partition.marginal_data().min_branch_length;
  let pass = GraphPass::new(graph)?;
  let nodes = &partition.marginal_data().nodes;
  let edges = &partition.marginal_data().edges;
  let outputs = pass.map_forward(
    nodes,
    edges,
    |key| treetime_utils::make_internal_error!("Partition node {key} is missing before the marginal forward pass"),
    |context| marginal_process_node_forward_indexed(partition, &gtr, min_branch_length, branch_lengths, &context),
  )?;
  partition.marginal_data_mut().nodes = outputs.nodes;
  partition.marginal_data_mut().edges = outputs.edges;
  Ok(())
}

fn marginal_process_node_forward_indexed(
  partition: &dyn IndexedMarginalPartition,
  gtr: &GTR,
  min_branch_length: f64,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  context: &GraphPassForwardContext<'_, DenseNodePartition, DenseEdgePartition, DenseNodePartition>,
) -> Result<GraphPassNodeOutput<DenseNodePartition, DenseEdgePartition>, Report> {
  let mut node = context.input.clone();

  // Clone this node's parent-edge input and overwrite only `msg_to_child`, exactly as before. The
  // edge's other fields (indels, msg_from_child, msg_to_parent, transmission) carry backward-pass state
  // and must survive the forward pass unchanged.
  let mut parent_edge = context.parent_edge.map(|(edge_key, edge)| (edge_key, edge.clone()));
  if let Some((edge_key, edge)) = parent_edge.as_mut() {
    let parent = context.parent.expect("Non-root node must have a parent");
    let safe_child = edge.msg_from_child.dis.mapv(|value| value.max(f64::MIN_POSITIVE));
    let mut dis = &parent.profile.dis / &safe_child;
    let delta_ll = normalize_inplace(&mut dis);
    let log_lh = forward_log_lh_remove_child(parent.profile.log_lh, edge.msg_from_child.log_lh);
    let log_lh = forward_log_lh_add_normalization(log_lh, delta_ll);
    edge.msg_to_child = DenseSeqDistribution { dis, log_lh };
    let branch_length = branch_lengths[&*edge_key].max(min_branch_length);
    let msg_child = gtr.evolve(&edge.msg_to_child.dis, branch_length, false);
    let mut dis = &edge.msg_to_parent.dis * &msg_child;
    let delta_ll = normalize_inplace(&mut dis);
    node.profile = DenseSeqDistribution {
      dis,
      log_lh: edge.msg_to_parent.log_lh + edge.msg_to_child.log_lh + LogLh::new(delta_ll),
    };
  }

  partition.indexed_forward_post(
    context.is_root,
    context.is_leaf,
    context.parent,
    &mut node,
    parent_edge.as_mut().map(|(_, edge)| edge),
  )?;

  let parent_message = parent_edge.map(|(_, edge)| edge);
  Ok(GraphPassNodeOutput { node, parent_message })
}
