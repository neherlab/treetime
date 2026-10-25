use crate::partition::assembly::{AlignmentGroup, alignment_groups};
use crate::partition::marginal::reconstruction::SparseReconstruction;
use crate::partition::storage::sparse::SparseEdgeObs;
use crate::seq::indel::InDel;
use eyre::Report;
use itertools::{Itertools, izip};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::reroot::EdgeMergeInfo;

pub(crate) fn assembled_edge_indels(
  group: &AlignmentGroup,
  partitions: &[SparseReconstruction],
  edge_key: GraphEdgeKey,
) -> Result<Vec<InDel>, Report> {
  group.layout.assemble_indels(&edge_pieces(group, partitions, edge_key))
}

pub(crate) fn store_edge_indels(
  group: &AlignmentGroup,
  partitions: &mut [SparseReconstruction],
  edge_key: GraphEdgeKey,
  indels: &[InDel],
) -> Result<(), Report> {
  let pieces = group.layout.split_indels(indels)?;
  store_edge_pieces(group, partitions, edge_key, pieces);
  Ok(())
}

pub(crate) fn merge_sparse_edges(
  partitions: &mut [SparseReconstruction],
  merges: &[EdgeMergeInfo],
) -> Result<(), Report> {
  let groups = sparse_alignment_groups(partitions)?;
  for merge in merges {
    merge_edge_pair(partitions, &groups, merge)?;
  }
  Ok(())
}

pub(crate) fn sparse_alignment_groups(partitions: &[SparseReconstruction]) -> Result<Vec<AlignmentGroup>, Report> {
  alignment_groups(
    partitions
      .iter()
      .map(|partition| (partition.layout(), partition.partition.length)),
  )
}

pub(crate) fn merge_edge_pair(
  partitions: &mut [SparseReconstruction],
  groups: &[AlignmentGroup],
  merge: &EdgeMergeInfo,
) -> Result<(), Report> {
  let pieces: Vec<_> = groups
    .iter()
    .map(|group| composed_edge_pieces(group, partitions, merge.parent_edge_key, merge.child_edge_key))
    .try_collect()?;
  for SparseReconstruction { partition, .. } in partitions.iter_mut() {
    let parent_edge = partition.obs_edges.remove(merge.parent_edge_key).unwrap_or_default();
    let child_edge = partition.obs_edges.remove(merge.child_edge_key).unwrap_or_default();
    let mut merged_edge = SparseEdgeObs::default();
    merged_edge.set_fitch_subs(parent_edge.chain_fitch_subs(child_edge.fitch_subs())?);
    partition.obs_edges.insert(merge.merged_edge_key, merged_edge);
  }
  for (group, pieces) in izip!(groups, pieces) {
    store_edge_pieces(group, partitions, merge.merged_edge_key, pieces);
  }
  Ok(())
}

pub(crate) fn composed_edge_pieces(
  group: &AlignmentGroup,
  partitions: &[SparseReconstruction],
  parent_edge_key: GraphEdgeKey,
  child_edge_key: GraphEdgeKey,
) -> Result<Vec<Vec<InDel>>, Report> {
  group.layout.compose_indels(
    &edge_pieces(group, partitions, parent_edge_key),
    &edge_pieces(group, partitions, child_edge_key),
  )
}

pub(crate) fn store_edge_pieces(
  group: &AlignmentGroup,
  partitions: &mut [SparseReconstruction],
  edge_key: GraphEdgeKey,
  pieces: Vec<Vec<InDel>>,
) {
  for (&index, pieces) in izip!(&group.partitions, pieces) {
    partitions[index]
      .partition
      .obs_edges
      .entry(edge_key)
      .or_default()
      .indels = pieces;
  }
}

fn edge_pieces(group: &AlignmentGroup, partitions: &[SparseReconstruction], edge_key: GraphEdgeKey) -> Vec<Vec<InDel>> {
  group
    .members(partitions)
    .into_iter()
    .map(|partition| {
      partition
        .partition
        .obs_edges
        .get(edge_key)
        .map_or_else(Vec::new, |edge| edge.indels.clone())
    })
    .collect()
}
