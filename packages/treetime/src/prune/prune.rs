use crate::optimize::topology::collapse::collapse_edge;
use crate::partition::marginal_dense::PartitionMarginalDense;
use crate::partition::marginal_sparse::PartitionMarginalSparse;
use crate::payload::ancestral::GraphAncestral;
use eyre::Report;
use itertools::Itertools;
use log::debug;
use parking_lot::RwLock;
use std::collections::BTreeSet;
use std::sync::Arc;
use treetime_graph::edge::{GraphEdgeKey, HasBranchLength};
use treetime_graph::node::{GraphNodeKey, Named};

pub fn prune_nodes(
  graph: &mut GraphAncestral,
  partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
  prune_short: Option<f64>,
  prune_empty: bool,
  node_names: &BTreeSet<String>,
) -> Result<(), Report> {
  prune_internal_nodes(graph, partitions, prune_short, prune_empty, node_names)?;
  graph.build()?;
  prune_leaves(graph, partitions, node_names)?;
  graph.build()?;
  Ok(())
}

/// Count current nucleotide mutations on one edge across all partitions.
pub fn get_edge_num_muts(
  partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
  edge_key: GraphEdgeKey,
) -> Result<Option<usize>, Report> {
  let mut total_muts = 0;
  let mut found_any = false;

  for partition in partitions {
    let partition = partition.read_arc();
    if let Some(edge) = partition.edges.get(&edge_key) {
      total_muts += edge.fitch_subs().len();
      found_any = true;
    }
  }

  Ok(found_any.then_some(total_muts))
}

pub fn collapse_sparse_edges_from_leaf_recursive(
  graph: &mut GraphAncestral,
  partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
  edge_key: GraphEdgeKey,
) -> Result<(), Report> {
  let mut current_edge_key = edge_key;
  let no_dense: &[Arc<RwLock<PartitionMarginalDense>>] = &[];

  loop {
    let parent_node_key = graph.get_source_node_key(current_edge_key)?;

    collapse_edge(graph, partitions, no_dense, current_edge_key)?;

    let next_edge_key = if should_collapse_parent(graph, parent_node_key) {
      graph.parent_inbound_edge(parent_node_key)?
    } else {
      None
    };

    match next_edge_key {
      Some(key) => current_edge_key = key,
      None => break,
    }
  }

  Ok(())
}

fn prune_internal_nodes(
  graph: &mut GraphAncestral,
  partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
  prune_short: Option<f64>,
  prune_empty: bool,
  node_names: &BTreeSet<String>,
) -> Result<(), Report> {
  #[allow(clippy::needless_collect)]
  let edges_to_collapse: Vec<_> = graph
    .get_edges()
    .iter()
    .map(|edge| -> Result<Option<GraphEdgeKey>, Report> {
      let edge = edge.read_arc();
      let target_is_leaf = graph.is_leaf(edge.target());

      if target_is_leaf {
        return Ok(None);
      }

      let weight = edge.payload().read_arc().branch_length();
      let should_prune_short = matches!((prune_short, weight), (Some(threshold), Some(weight)) if weight < threshold);

      let should_prune_empty = prune_empty && get_edge_num_muts(partitions, edge.key())? == Some(0);

      let target_node = graph
        .get_node(edge.target())
        .map(|n| n.read_arc().payload().read_arc().name().map(|n| n.as_ref().to_owned()));
      let should_prune_by_name = target_node.is_some_and(|name| name.is_some_and(|n| node_names.contains(&n)));

      let should_prune = should_prune_short || should_prune_empty || should_prune_by_name;
      Ok(should_prune.then(|| edge.key()))
    })
    .collect::<Result<Vec<_>, Report>>()?
    .into_iter()
    .flatten()
    .collect();

  let no_dense: &[Arc<RwLock<PartitionMarginalDense>>] = &[];
  edges_to_collapse.into_iter().try_for_each(|edge_key| {
    debug!("Collapsing internal edge: {edge_key}");
    collapse_edge(graph, partitions, no_dense, edge_key)
  })
}

fn prune_leaves(
  graph: &mut GraphAncestral,
  partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
  node_names: &BTreeSet<String>,
) -> Result<(), Report> {
  #[allow(clippy::needless_collect)]
  let edges_to_collapse = graph
    .get_edges()
    .iter()
    .filter_map(|edge| {
      let edge = edge.read_arc();
      let target_is_leaf = graph.is_leaf(edge.target());

      if !target_is_leaf {
        return None;
      }

      let target_node = graph.get_node(edge.target())?.read_arc().payload().read_arc();
      let name = target_node.name();
      let should_prune_by_name = name.is_some_and(|name| node_names.contains(name.as_ref()));

      should_prune_by_name.then(|| edge.key())
    })
    .collect_vec();

  edges_to_collapse.into_iter().try_for_each(|edge_key| {
    debug!("Collapsing leaf edge: {edge_key}");
    collapse_sparse_edges_from_leaf_recursive(graph, partitions, edge_key)
  })
}

fn should_collapse_parent(graph: &GraphAncestral, node_key: GraphNodeKey) -> bool {
  graph.has_at_most_one_child(node_key) && !graph.is_root(node_key)
}
