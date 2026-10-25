use crate::gtr::jc_distance::jukes_cantor_distance;
use crate::optimize::topology::polytomy_nodes::find_polytomy_nodes;
use crate::partition::assembly::AlignmentGroup;
use crate::partition::marginal::reconstruction::SparseReconstruction;
use crate::partition::marginal::sparse::alignment_indels::{
  assembled_edge_indels, sparse_alignment_groups, store_edge_indels,
};
use crate::partition::storage::sparse::SparseNodeObs;
use crate::seq::indel::InDel;
use crate::seq::mutation::Sub;
use eyre::Report;
use itertools::{Itertools, izip};
use log::debug;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::EdgeMap;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::node_names::NodeNames;
use treetime_utils::cast::ToF64Exact;
use treetime_utils::iterator::difference::iterator_difference;

pub(crate) fn merge_shared_mutation_branches(
  graph: &mut Graph,
  names: &mut NodeNames,
  partitions: &mut [SparseReconstruction],
  branch_lengths: &mut EdgeMap<Option<f64>>,
) -> Result<usize, Report> {
  let alignments = sparse_alignment_groups(partitions)?;
  let mut total_merged = 0;

  let polytomy_keys = find_polytomy_nodes(graph);

  for node_key in polytomy_keys {
    let merged = merge_single_polytomy(graph, names, partitions, &alignments, node_key, branch_lengths)?;
    total_merged += merged;
  }

  if total_merged > 0 {
    debug!("Merged {total_merged} sibling groups sharing mutations");
  } else {
    debug!("No sibling branches sharing mutations found");
  }

  Ok(total_merged)
}

pub(crate) fn merge_single_polytomy(
  graph: &mut Graph,
  names: &mut NodeNames,
  partitions: &mut [SparseReconstruction],
  alignments: &[AlignmentGroup],
  node_key: GraphNodeKey,
  branch_lengths: &mut EdgeMap<Option<f64>>,
) -> Result<usize, Report> {
  let mut nodes_created = 0;

  loop {
    let child_edges = collect_child_edge_keys(graph, node_key);
    if child_edges.len() <= 2 {
      break;
    }

    let index = build_mutation_index(partitions, alignments, &child_edges)?;
    let candidates = find_merge_groups(&index, partitions.len(), alignments.len());
    let selected = greedy_disjoint_group_matching(candidates);

    if selected.is_empty() {
      break;
    }

    for group in selected {
      let edges_fmt = group.edges.iter().format(", ");
      let shared_subs_fmt = group.shared_subs.iter().flatten().format(", ");
      let shared_indels_fmt = group.shared_indels.iter().flatten().format(", ");
      debug!(
        "Merging {} siblings under node {node_key}: edges [{edges_fmt}] share {} mutations: subs=[{shared_subs_fmt}] indels=[{shared_indels_fmt}]",
        group.edges.len(),
        group.total_shared,
      );
      merge_sibling_group(graph, names, partitions, alignments, node_key, &group, branch_lengths)?;
      nodes_created += 1;
    }
  }

  Ok(nodes_created)
}

fn collect_child_edge_keys(graph: &Graph, node_key: GraphNodeKey) -> Vec<GraphEdgeKey> {
  let node = &graph[node_key];
  node.outbound().to_vec()
}

fn build_mutation_index(
  partitions: &[SparseReconstruction],
  alignments: &[AlignmentGroup],
  child_edges: &[GraphEdgeKey],
) -> Result<MutationIndex, Report> {
  let mut index: MutationIndex = BTreeMap::new();
  for &edge_key in child_edges {
    for (pi, SparseReconstruction { partition, .. }) in partitions.iter().enumerate() {
      let empty_subs: &[Sub] = &[];
      let edge = partition.obs_edges.get(edge_key);
      for sub in edge.map_or(empty_subs, |e| e.fitch_subs()) {
        index
          .entry(MutationKey::Sub(pi, sub.clone()))
          .or_default()
          .push(edge_key);
      }
    }
    for (ai, alignment) in alignments.iter().enumerate() {
      for indel in assembled_edge_indels(alignment, partitions, edge_key)? {
        index.entry(MutationKey::InDel(ai, indel)).or_default().push(edge_key);
      }
    }
  }
  Ok(index)
}

fn find_merge_groups(index: &MutationIndex, n_partitions: usize, n_alignments: usize) -> Vec<MergeGroup> {
  let mut seen: BTreeSet<Vec<GraphEdgeKey>> = BTreeSet::new();
  let mut groups = Vec::new();

  for edges in index.values() {
    if edges.len() < 2 {
      continue;
    }
    let mut group_edges = edges.clone();
    group_edges.sort_unstable();
    group_edges.dedup();
    if !seen.insert(group_edges.clone()) {
      continue;
    }
    let (shared_subs, shared_indels) = shared_mutations_for_group(index, n_partitions, n_alignments, &group_edges);
    let total_shared: usize =
      shared_subs.iter().map(|v| v.len()).sum::<usize>() + shared_indels.iter().map(|v| v.len()).sum::<usize>();
    if total_shared > 0 {
      groups.push(MergeGroup {
        edges: group_edges,
        shared_subs,
        shared_indels,
        total_shared,
      });
    }
  }
  groups
}

fn shared_mutations_for_group(
  index: &MutationIndex,
  n_partitions: usize,
  n_alignments: usize,
  group: &[GraphEdgeKey],
) -> (Vec<Vec<Sub>>, Vec<Vec<InDel>>) {
  let mut shared_subs: Vec<Vec<Sub>> = vec![Vec::new(); n_partitions];
  let mut shared_indels: Vec<Vec<InDel>> = vec![Vec::new(); n_alignments];
  for (key, edges) in index {
    if group.iter().all(|e| edges.contains(e)) {
      match key {
        MutationKey::Sub(pi, sub) => shared_subs[*pi].push(sub.clone()),
        MutationKey::InDel(ai, indel) => shared_indels[*ai].push(indel.clone()),
      }
    }
  }
  (shared_subs, shared_indels)
}

type MutationIndex = BTreeMap<MutationKey, Vec<GraphEdgeKey>>;

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
enum MutationKey {
  Sub(usize, Sub),
  InDel(usize, InDel),
}

fn greedy_disjoint_group_matching(mut groups: Vec<MergeGroup>) -> Vec<MergeGroup> {
  groups.sort_unstable_by(|a, b| {
    b.total_shared
      .cmp(&a.total_shared)
      .then(b.edges.len().cmp(&a.edges.len()))
      .then_with(|| a.edges.cmp(&b.edges))
  });
  let mut used: BTreeSet<GraphEdgeKey> = BTreeSet::new();
  groups.retain(|g| {
    if g.edges.iter().all(|e| !used.contains(e)) {
      used.extend(g.edges.iter().copied());
      true
    } else {
      false
    }
  });
  groups
}

fn merge_sibling_group(
  graph: &mut Graph,
  names: &mut NodeNames,
  partitions: &mut [SparseReconstruction],
  alignments: &[AlignmentGroup],
  parent_key: GraphNodeKey,
  group: &MergeGroup,
  branch_lengths: &mut EdgeMap<Option<f64>>,
) -> Result<(), Report> {
  let child_keys: Vec<GraphNodeKey> = group
    .edges
    .iter()
    .map(|&ek| graph.get_target_node_key(ek))
    .collect::<Result<Vec<_>, _>>()?;

  let Remaining {
    subs: remaining_subs,
    indels: remaining_indels,
  } = remaining_mutations(partitions, alignments, group)?;

  let total_alignment_length: usize = partitions.iter().map(|p| p.partition.length).sum();
  let n_states = partitions[0].partition.alphabet.n_canonical();
  let jc_bl = |count: usize| -> f64 {
    if total_alignment_length > 0 {
      jukes_cantor_distance(count.to_f64_exact() / total_alignment_length.to_f64_exact(), n_states)
    } else {
      0.0
    }
  };

  let new_edge_bl = jc_bl(group.total_shared);
  let child_bls: Vec<f64> = (0..group.edges.len())
    .map(|i| {
      let remaining = remaining_subs.iter().map(|subs| subs[i].len()).sum::<usize>()
        + remaining_indels.iter().map(|indels| indels[i].len()).sum::<usize>();
      jc_bl(remaining)
    })
    .collect();

  for &ek in &group.edges {
    graph.remove_edge(ek)?;
    branch_lengths.remove(ek);
  }

  let new_node_key = graph.add_node();
  names.name_new_node(new_node_key);

  let new_parent_edge_key = graph.add_edge(parent_key, new_node_key)?;
  branch_lengths.insert(new_parent_edge_key, Some(new_edge_bl));

  let mut new_child_edge_keys: Vec<GraphEdgeKey> = Vec::with_capacity(child_keys.len());
  for (&ck, &bl) in child_keys.iter().zip(child_bls.iter()) {
    let new_ek = graph.add_edge(new_node_key, ck)?;
    branch_lengths.insert(new_ek, Some(bl));
    new_child_edge_keys.push(new_ek);
  }

  for (pi, SparseReconstruction { partition, .. }) in partitions.iter_mut().enumerate() {
    let mut new_node = SparseNodeObs::empty(&partition.alphabet);
    new_node.composition = partition.obs_nodes[&parent_key].composition.clone();
    partition.obs_nodes.entry(new_node_key).or_insert(new_node);

    for &ek in &group.edges {
      partition.obs_edges.remove(ek);
    }

    let parent_edge = partition.obs_edges.entry(new_parent_edge_key).or_default();
    parent_edge.set_fitch_subs(group.shared_subs[pi].clone());

    for (i, &new_ek) in new_child_edge_keys.iter().enumerate() {
      let child_edge = partition.obs_edges.entry(new_ek).or_default();
      child_edge.set_fitch_subs(remaining_subs[pi][i].clone());
    }
  }
  for (alignment, shared_indels, remaining_indels) in izip!(alignments, &group.shared_indels, &remaining_indels) {
    store_edge_indels(alignment, partitions, new_parent_edge_key, shared_indels)?;
    for (&new_ek, remaining_indels) in izip!(&new_child_edge_keys, remaining_indels) {
      store_edge_indels(alignment, partitions, new_ek, remaining_indels)?;
    }
  }

  Ok(())
}

fn remaining_mutations(
  partitions: &[SparseReconstruction],
  alignments: &[AlignmentGroup],
  group: &MergeGroup,
) -> Result<Remaining, Report> {
  let remaining_subs: Vec<Vec<Vec<Sub>>> = partitions
    .iter()
    .enumerate()
    .map(|(pi, SparseReconstruction { partition, .. })| {
      let empty_subs: &[Sub] = &[];
      group
        .edges
        .iter()
        .map(|&ek| {
          let all_subs = partition.obs_edges.get(ek).map_or(empty_subs, |e| e.fitch_subs());
          iterator_difference(all_subs, &group.shared_subs[pi])
            .cloned()
            .collect_vec()
        })
        .collect()
    })
    .collect();
  let remaining_indels: Vec<Vec<Vec<InDel>>> = izip!(alignments, &group.shared_indels)
    .map(|(alignment, shared_indels)| {
      group
        .edges
        .iter()
        .map(|&ek| {
          let all_indels = assembled_edge_indels(alignment, partitions, ek)?;
          Ok(iterator_difference(all_indels, shared_indels.iter().cloned()).collect_vec())
        })
        .collect::<Result<Vec<_>, Report>>()
    })
    .collect::<Result<Vec<_>, Report>>()?;
  Ok(Remaining {
    subs: remaining_subs,
    indels: remaining_indels,
  })
}

struct MergeGroup {
  edges: Vec<GraphEdgeKey>,
  shared_subs: Vec<Vec<Sub>>,
  shared_indels: Vec<Vec<InDel>>,
  total_shared: usize,
}

struct Remaining {
  subs: Vec<Vec<Vec<Sub>>>,
  indels: Vec<Vec<Vec<InDel>>>,
}
