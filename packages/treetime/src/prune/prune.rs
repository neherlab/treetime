use crate::optimize::topology::collapse::collapse_edge;
use crate::partition::assembly::AlignmentGroup;
use crate::partition::marginal::reconstruction::SparseReconstruction;
use crate::partition::marginal::sparse::alignment_indels::sparse_alignment_groups;
use eyre::Report;
use itertools::Itertools;
use log::debug;
use std::collections::BTreeSet;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::EdgeMap;
use treetime_graph::node::GraphNodeKey;

pub(crate) fn prune_nodes(
  graph: &mut Graph,
  partitions: &mut [SparseReconstruction],
  prune_short: Option<f64>,
  prune_empty: bool,
  node_keys: &BTreeSet<GraphNodeKey>,
  branch_lengths: &mut EdgeMap<Option<f64>>,
) -> Result<(), Report> {
  let groups = sparse_alignment_groups(partitions)?;
  prune_internal_nodes(
    graph,
    partitions,
    &groups,
    prune_short,
    prune_empty,
    node_keys,
    branch_lengths,
  )?;
  graph.build()?;
  prune_leaves(graph, partitions, &groups, node_keys, branch_lengths)?;
  graph.build()?;
  Ok(())
}

fn prune_internal_nodes(
  graph: &mut Graph,
  partitions: &mut [SparseReconstruction],
  groups: &[AlignmentGroup],
  prune_short: Option<f64>,
  prune_empty: bool,
  node_keys: &BTreeSet<GraphNodeKey>,
  branch_lengths: &mut EdgeMap<Option<f64>>,
) -> Result<(), Report> {
  #[expect(
    clippy::needless_collect,
    reason = "collecting ends the graph borrow before edges are collapsed"
  )]
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

      let should_prune_by_key = node_keys.contains(&edge.target());

      let should_prune = should_prune_short || should_prune_empty || should_prune_by_key;
      Ok(should_prune.then(|| edge.key()))
    })
    .collect::<Result<Vec<_>, Report>>()?
    .into_iter()
    .flatten()
    .collect();

  edges_to_collapse.into_iter().try_for_each(|edge_key| {
    debug!("Collapsing internal edge: {edge_key}");
    collapse_edge(graph, partitions, groups, edge_key, branch_lengths)
  })
}

pub(crate) fn get_edge_num_muts(
  partitions: &[SparseReconstruction],
  edge_key: GraphEdgeKey,
) -> Result<Option<usize>, Report> {
  let mut total_muts = 0;
  let mut found_any = false;

  for reconstruction in partitions {
    if let Some(edge) = reconstruction.partition.obs_edges.get(edge_key) {
      total_muts += edge.fitch_subs().len();
      found_any = true;
    }
  }

  Ok(found_any.then_some(total_muts))
}

fn prune_leaves(
  graph: &mut Graph,
  partitions: &mut [SparseReconstruction],
  groups: &[AlignmentGroup],
  node_keys: &BTreeSet<GraphNodeKey>,
  branch_lengths: &mut EdgeMap<Option<f64>>,
) -> Result<(), Report> {
  let edges_to_collapse = graph
    .get_edges()
    .filter_map(|edge| {
      let target_is_leaf = graph.is_leaf(edge.target());

      if !target_is_leaf {
        return None;
      }

      node_keys.contains(&edge.target()).then(|| edge.key())
    })
    .collect_vec();

  edges_to_collapse.into_iter().try_for_each(|edge_key| {
    debug!("Collapsing leaf edge: {edge_key}");
    collapse_sparse_edges_from_leaf_recursive(graph, partitions, groups, edge_key, branch_lengths)
  })
}

pub(crate) fn collapse_sparse_edges_from_leaf_recursive(
  graph: &mut Graph,
  partitions: &mut [SparseReconstruction],
  groups: &[AlignmentGroup],
  edge_key: GraphEdgeKey,
  branch_lengths: &mut EdgeMap<Option<f64>>,
) -> Result<(), Report> {
  let mut current_edge_key = edge_key;

  loop {
    let parent_node_key = graph.get_source_node_key(current_edge_key)?;

    collapse_edge(graph, partitions, groups, current_edge_key, branch_lengths)?;

    let next_edge_key = if should_collapse_parent(graph, parent_node_key) {
      graph.node_parent(parent_node_key)?.map(|(_, edge_key)| edge_key)
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

fn should_collapse_parent(graph: &Graph, node_key: GraphNodeKey) -> bool {
  graph.has_at_most_one_child(node_key) && !graph.is_root(node_key)
}
