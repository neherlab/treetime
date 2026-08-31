use crate::partition::marginal::shared::data::IndexedMarginalPartition;
use crate::partition::marginal::shared::normalize::{
  forward_log_lh_add_normalization, forward_log_lh_remove_child, normalize_from_log, normalize_inplace,
};
use crate::partition::storage::dense::{DenseEdgePartition, DenseNodePartition, DenseSeqDistribution};
use eyre::Report;
use itertools::Itertools;
use std::collections::BTreeMap;
use treetime_graph::edge::EdgeOptimizeOps;
use treetime_graph::graph::Graph;
use treetime_graph::indexed_pass::{IndexedPass, IndexedPassDependencies, IndexedPassSlot};
use treetime_graph::node::{GraphNode, Named};
use treetime_primitives::LogLh;

pub fn marginal_process_backward_indexed<N, E>(
  partition: &mut impl IndexedMarginalPartition<N, E>,
  graph: &Graph<N, E, ()>,
) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  let mut missing_nodes = BTreeMap::new();
  for node in graph.get_nodes() {
    let key = node.read_arc().key();
    if !partition.marginal_data().nodes.contains_key(&key) {
      missing_nodes.insert(key, partition.indexed_missing_node(key)?);
    }
  }
  partition.marginal_data_mut().nodes.append(&mut missing_nodes);
  let (nodes, edges) = partition.indexed_storage_mut();
  let mut pass = IndexedPass::new(graph, nodes, edges, |_| unreachable!("Missing nodes were initialized"))?;
  let result = pass.try_for_each_backward(|dependencies, slot| {
    marginal_process_node_backward_indexed(partition, graph, dependencies, slot)
  });
  let (nodes, edges) = pass.into_maps()?;
  partition.marginal_data_mut().nodes = nodes;
  partition.marginal_data_mut().edges = edges;
  result
}

fn marginal_process_node_backward_indexed<N, E>(
  partition: &impl IndexedMarginalPartition<N, E>,
  graph: &Graph<N, E, ()>,
  dependencies: &IndexedPassDependencies<DenseNodePartition, DenseEdgePartition>,
  slot: &mut IndexedPassSlot<DenseNodePartition, DenseEdgePartition>,
) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  let graph_node = graph.get_node(slot.key).expect("Indexed node must exist in graph");
  let graph_node = graph_node.read_arc();
  let msg_to_parent = if graph_node.is_leaf() {
    partition.indexed_leaf_profile(&slot.node)?
  } else {
    let child_pairs = graph.children_of(&graph_node);
    let children = child_pairs
      .iter()
      .map(|(child, _)| {
        let child_key = child.read_arc().key();
        dependencies.node(child_key)
      })
      .collect_vec();
    slot.node.seq = partition.indexed_backward_internal(&children)?;

    let child_edges = child_pairs
      .iter()
      .map(|(_, edge)| {
        let edge_key = edge.read_arc().key();
        dependencies.edge(edge_key)
      })
      .collect_vec();
    let first_edge = child_edges.first().expect("Internal node must have children");
    let mut log_dis = first_edge.msg_from_child.dis.mapv(f64::ln);
    for edge in child_edges.iter().skip(1) {
      log_dis += &edge.msg_from_child.dis.mapv(f64::ln);
    }
    let (dis, delta_ll) = normalize_from_log(&log_dis);
    let log_lh = child_edges.iter().map(|edge| edge.msg_from_child.log_lh).sum::<LogLh>() + LogLh::new(delta_ll);
    slot.node.profile = DenseSeqDistribution {
      dis: dis.clone(),
      log_lh,
    };
    DenseSeqDistribution { dis, log_lh }
  };

  if graph_node.is_root() {
    let data = partition.marginal_data();
    let mut dis = &msg_to_parent.dis * &data.gtr.pi;
    let delta_ll = normalize_inplace(&mut dis);
    slot.node.profile = DenseSeqDistribution {
      dis,
      log_lh: msg_to_parent.log_lh + LogLh::new(delta_ll),
    };
  } else {
    let (edge_key, edge) = slot
      .parent_edge
      .as_mut()
      .expect("Non-root node must own its parent edge");
    let branch_length = graph
      .get_edge(*edge_key)
      .expect("Indexed edge must exist in graph")
      .read_arc()
      .payload()
      .read_arc()
      .profile_branch_length()
      .unwrap_or(0.0);
    let branch_length = partition.marginal_data().effective_branch_length(branch_length);
    edge.msg_from_child = DenseSeqDistribution {
      dis: partition
        .marginal_data()
        .gtr
        .propagate_profile(&msg_to_parent.dis, branch_length, false),
      log_lh: msg_to_parent.log_lh,
    };
    edge.msg_to_parent = msg_to_parent;
  }
  Ok(())
}

pub fn marginal_process_forward_indexed<N, E>(
  partition: &mut impl IndexedMarginalPartition<N, E>,
  graph: &Graph<N, E, ()>,
) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  let (nodes, edges) = partition.indexed_storage_mut();
  let mut pass = IndexedPass::new(graph, nodes, edges, |key| {
    treetime_utils::make_internal_error!("Partition node {key} is missing before the marginal forward pass")
  })?;
  let result = pass.try_for_each_forward(|dependencies, slot| {
    marginal_process_node_forward_indexed(partition, graph, dependencies, slot)
  });
  let (nodes, edges) = pass.into_maps()?;
  partition.marginal_data_mut().nodes = nodes;
  partition.marginal_data_mut().edges = edges;
  result
}

#[allow(clippy::too_many_arguments)]
fn marginal_process_node_forward_indexed<N, E>(
  partition: &impl IndexedMarginalPartition<N, E>,
  graph: &Graph<N, E, ()>,
  dependencies: &IndexedPassDependencies<DenseNodePartition, DenseEdgePartition>,
  slot: &mut IndexedPassSlot<DenseNodePartition, DenseEdgePartition>,
) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  let mut parent_edge = None;
  let parent = if let Some(parent_key) = slot.parent_key {
    let parent = dependencies.node(parent_key);
    parent_edge = slot.parent_edge.as_mut();
    Some(parent)
  } else {
    None
  };

  if let Some((edge_key, edge)) = parent_edge.as_mut() {
    let safe_child = edge.msg_from_child.dis.mapv(|value| value.max(f64::MIN_POSITIVE));
    let mut dis = &parent.expect("Non-root node must have a parent").profile.dis / &safe_child;
    let delta_ll = normalize_inplace(&mut dis);
    let log_lh = forward_log_lh_remove_child(
      parent.expect("Non-root node must have a parent").profile.log_lh,
      edge.msg_from_child.log_lh,
    );
    let log_lh = forward_log_lh_add_normalization(log_lh, delta_ll);
    edge.msg_to_child = DenseSeqDistribution { dis, log_lh };
    let branch_length = graph
      .get_edge(*edge_key)
      .expect("Indexed edge must exist in graph")
      .read_arc()
      .payload()
      .read_arc()
      .profile_branch_length()
      .unwrap_or(0.0);
    let branch_length = partition.marginal_data().effective_branch_length(branch_length);
    let msg_child = partition
      .marginal_data()
      .gtr
      .evolve(&edge.msg_to_child.dis, branch_length, false);
    let mut dis = &edge.msg_to_parent.dis * &msg_child;
    let delta_ll = normalize_inplace(&mut dis);
    slot.node.profile = DenseSeqDistribution {
      dis,
      log_lh: edge.msg_to_parent.log_lh + edge.msg_to_child.log_lh + LogLh::new(delta_ll),
    };
  }

  partition.indexed_forward_post(
    slot.parent_key.is_none(),
    graph.is_leaf(slot.key),
    parent,
    &mut slot.node,
    parent_edge.as_mut().map(|(_, edge)| edge),
  )?;
  Ok(())
}
