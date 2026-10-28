use crate::branch_lengths::one_mutation;
use crate::optimize::indel::EdgeIndels;
use crate::partition::assembly::AlignmentGroup;
use crate::partition::marginal::reconstruction::SparseReconstruction;
use crate::partition::marginal::sparse::alignment_indels::{
  composed_edge_pieces, sparse_alignment_groups, store_edge_pieces,
};
use crate::partition::set::PartitionSet;
use crate::seq::mutation::compose_substitutions;
use eyre::Report;
use itertools::Itertools;
use log::debug;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::EdgeMap;
use treetime_graph::schedule::GraphSchedule;
use treetime_utils::interval::range::range_contains;

const SHORT_BRANCH_MUTATIONS: f64 = 0.1;

pub(crate) fn collapse_edge(
  graph: &mut Graph,
  sparse: &mut [SparseReconstruction],
  groups: &[AlignmentGroup],
  edge_key: GraphEdgeKey,
  branch_lengths: &mut EdgeMap<Option<f64>>,
) -> Result<(), Report> {
  let target_node_key = graph.get_target_node_key(edge_key)?;

  let removed_bl = branch_lengths[&edge_key];
  let (_, _removed_edge, new_edges) = graph.collapse_edge(edge_key)?;

  for new_edge in &new_edges {
    let new_edge_key = *new_edge;

    if let (Some(bl1), Some(bl2)) = (removed_bl, branch_lengths[&new_edge_key]) {
      branch_lengths.insert(new_edge_key, Some(bl1 + bl2));
    }

    let child_key = graph.get_target_node_key(new_edge_key)?;
    for reconstruction in sparse.iter_mut() {
      let partition = &mut reconstruction.partition;
      let child_non_char = partition
        .obs_nodes
        .get(child_key)
        .map_or_else(Vec::new, |child| child.non_char.clone());
      let removed_edge_data = partition.obs_edges[&edge_key].clone();
      let child_edge = partition.obs_edges.entry(new_edge_key).or_default();
      let merged_subs = removed_edge_data
        .chain_fitch_subs(child_edge.fitch_subs())?
        .into_iter()
        .filter(|sub| !range_contains(&child_non_char, sub.pos()))
        .collect();
      child_edge.set_fitch_subs(merged_subs);

      let estimates = &mut reconstruction.edges.estimates;
      match (estimates.get(edge_key), estimates.get(new_edge_key)) {
        (Some(removed), Some(child)) => {
          let composed = compose_substitutions(removed, child)?;
          estimates.insert(new_edge_key, composed);
        },
        _ => {
          estimates.remove(new_edge_key);
        },
      }
    }
    for group in groups {
      let pieces = composed_edge_pieces(group, sparse, edge_key, new_edge_key)?;
      store_edge_pieces(group, sparse, new_edge_key, pieces);
    }
  }

  for reconstruction in sparse.iter_mut() {
    reconstruction.partition.obs_nodes.remove(target_node_key);
    reconstruction.partition.obs_edges.remove(edge_key);
    reconstruction.edges.estimates.remove(edge_key);
  }

  branch_lengths.remove(edge_key);

  Ok(())
}

pub(crate) fn collapse_short_internal_edges(
  graph: &mut Graph,
  mut partitions: PartitionSet,
  indels: &EdgeIndels,
  branch_lengths: &mut EdgeMap<Option<f64>>,
) -> Result<PartitionSet, Report> {
  let total_length = partitions.total_length();
  let threshold = if total_length > 0 {
    SHORT_BRANCH_MUTATIONS * one_mutation(total_length)
  } else {
    0.0
  };
  let schedule = GraphSchedule::new(graph)?;
  let short_edges: Vec<GraphEdgeKey> = schedule
    .postorder()
    .iter()
    .filter(|&&node_key| !graph.is_leaf(node_key))
    .flat_map(|&node_key| schedule.parents(node_key).iter().map(|&(_, edge_key)| edge_key))
    .filter(|&edge_key| !indels.term(edge_key).requires_positive_length())
    .map(|edge_key| -> Result<_, Report> {
      let length = branch_lengths[&edge_key].unwrap_or(0.0);
      let is_short = length == 0.0 || (length < threshold && !carries_substitution(graph, &partitions, edge_key)?);
      Ok(is_short.then_some(edge_key))
    })
    .flatten_ok()
    .try_collect()?;
  if short_edges.is_empty() {
    return Ok(partitions);
  }

  let sparse = partitions.sparse_mut();
  let groups = sparse_alignment_groups(sparse)?;
  for &edge_key in &short_edges {
    collapse_edge(graph, sparse, &groups, edge_key, branch_lengths)?;
  }
  debug!(
    "Collapsed {} internal edges shorter than {threshold} without substitutions after the optimization loop",
    short_edges.len()
  );
  graph.build()?;
  Ok(partitions.reconcile_node_states(graph))
}

fn carries_substitution(graph: &Graph, partitions: &PartitionSet, edge_key: GraphEdgeKey) -> Result<bool, Report> {
  partitions
    .iter()
    .map(|partition| partition.edge_sub_count(graph, edge_key))
    .process_results(|mut counts| counts.any(|count| count.is_none_or(|count| count > 0)))
}
