use crate::optimize::topology::collapse::collapse_edge;
use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
use eyre::Report;
use itertools::Itertools;
use log::debug;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;

pub fn prune_nodes(
  graph: &mut Graph,
  partitions: &mut [PartitionMarginalSparse],
  prune_short: Option<f64>,
  prune_empty: bool,
  node_names: &BTreeSet<String>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<(), Report> {
  // `names` is the pre-prune node-name map propagated from the command entry: collapse removes nodes
  // but never renames survivors, so the pre-prune label of every surviving node is its final label.
  // `branch_lengths` is the value map the collapse producers read and update in place, and it exits
  // reflecting the pruned tree.
  prune_internal_nodes(
    graph,
    partitions,
    prune_short,
    prune_empty,
    node_names,
    names,
    branch_lengths,
  )?;
  graph.build()?;
  prune_leaves(graph, partitions, node_names, names, branch_lengths)?;
  graph.build()?;
  Ok(())
}

/// Count current nucleotide mutations on one edge across all partitions.
pub fn get_edge_num_muts(
  partitions: &[PartitionMarginalSparse],
  edge_key: GraphEdgeKey,
) -> Result<Option<usize>, Report> {
  let mut total_muts = 0;
  let mut found_any = false;

  for partition in partitions {
    if let Some(edge) = partition.obs_edges.get(&edge_key) {
      total_muts += edge.fitch_subs().len();
      found_any = true;
    }
  }

  Ok(found_any.then_some(total_muts))
}

pub fn collapse_sparse_edges_from_leaf_recursive(
  graph: &mut Graph,
  partitions: &mut [PartitionMarginalSparse],
  edge_key: GraphEdgeKey,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<(), Report> {
  let mut current_edge_key = edge_key;

  loop {
    let parent_node_key = graph.get_source_node_key(current_edge_key)?;

    collapse_edge(graph, partitions, current_edge_key, branch_lengths)?;

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
  graph: &mut Graph,
  partitions: &mut [PartitionMarginalSparse],
  prune_short: Option<f64>,
  prune_empty: bool,
  node_names: &BTreeSet<String>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<(), Report> {
  #[allow(clippy::needless_collect)]
  let edges_to_collapse: Vec<_> = graph
    .get_edges()
    .map(|edge| -> Result<Option<GraphEdgeKey>, Report> {
      let target_is_leaf = graph.is_leaf(edge.target());

      if target_is_leaf {
        return Ok(None);
      }

      let weight = branch_lengths[&edge.key()];
      let should_prune_short = matches!((prune_short, weight), (Some(threshold), Some(weight)) if weight < threshold);

      let should_prune_empty = prune_empty && get_edge_num_muts(partitions, edge.key())? == Some(0);

      let should_prune_by_name = names[&edge.target()]
        .as_deref()
        .is_some_and(|name| node_names.contains(name));

      let should_prune = should_prune_short || should_prune_empty || should_prune_by_name;
      Ok(should_prune.then(|| edge.key()))
    })
    .collect::<Result<Vec<_>, Report>>()?
    .into_iter()
    .flatten()
    .collect();

  edges_to_collapse.into_iter().try_for_each(|edge_key| {
    debug!("Collapsing internal edge: {edge_key}");
    collapse_edge(graph, partitions, edge_key, branch_lengths)
  })
}

fn prune_leaves(
  graph: &mut Graph,
  partitions: &mut [PartitionMarginalSparse],
  node_names: &BTreeSet<String>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<(), Report> {
  #[allow(clippy::needless_collect)]
  let edges_to_collapse = graph
    .get_edges()
    .filter_map(|edge| {
      let target_is_leaf = graph.is_leaf(edge.target());

      if !target_is_leaf {
        return None;
      }

      let should_prune_by_name = names[&edge.target()]
        .as_deref()
        .is_some_and(|name| node_names.contains(name));

      should_prune_by_name.then(|| edge.key())
    })
    .collect_vec();

  edges_to_collapse.into_iter().try_for_each(|edge_key| {
    debug!("Collapsing leaf edge: {edge_key}");
    collapse_sparse_edges_from_leaf_recursive(graph, partitions, edge_key, branch_lengths)
  })
}

fn should_collapse_parent(graph: &Graph, node_key: GraphNodeKey) -> bool {
  graph.has_at_most_one_child(node_key) && !graph.is_root(node_key)
}
