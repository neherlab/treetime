use crate::gtr::jc_distance::jukes_cantor_distance;
use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
use crate::representation::payload::ancestral::{EdgeAncestral, GraphAncestral, NodeAncestral};
use crate::representation::payload::sparse::SparseNodePartition;
use crate::seq::indel::InDel;
use crate::seq::mutation::Sub;
use eyre::Report;
use itertools::Itertools;
use log::debug;
use parking_lot::RwLock;
use std::collections::{BTreeMap, HashSet};
use std::sync::Arc;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::node::GraphNodeKey;
use treetime_utils::iterator::difference::iterator_difference;

/// Merge sibling branches in polytomies that share identical substitutions.
///
/// Tree builders produce arbitrary binary resolutions of polytomies. When sibling branches
/// carry identical substitutions, they represent redundant resolutions. This function groups
/// such siblings (two or more) under a new internal node:
///
/// Before: P-A (subs: {G100T, A200C}), P-B (subs: {G100T, A200C, T300G}), P-C (subs: {G100T, A200C})
/// After:  P-N (subs: {G100T, A200C}), N-A (subs: {}), N-B (subs: {T300G}), N-C (subs: {})
///
/// Each round finds groups of k >= 2 siblings whose mutation sets share a non-empty
/// intersection, selects a disjoint set of groups, and merges each group under one new
/// internal node. Repeats until no siblings share any mutations. Returns the total number
/// of new internal nodes created.
pub fn merge_shared_mutation_branches(
  graph: &mut GraphAncestral,
  partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
) -> Result<usize, Report> {
  let mut total_merged = 0;

  let polytomy_keys = find_polytomy_nodes(graph);

  for node_key in polytomy_keys {
    let merged = merge_single_polytomy(graph, partitions, node_key)?;
    total_merged += merged;
  }

  if total_merged > 0 {
    debug!("Merged {total_merged} sibling groups sharing mutations");
  } else {
    debug!("No sibling branches sharing mutations found");
  }

  Ok(total_merged)
}

fn find_polytomy_nodes(graph: &GraphAncestral) -> Vec<GraphNodeKey> {
  graph
    .get_nodes()
    .into_iter()
    .filter_map(|node| {
      let node = node.read_arc();
      (node.degree_out() > 2).then_some(node.key())
    })
    .collect_vec()
}

/// Resolve one polytomy by greedy batch merging of siblings sharing mutations.
///
/// Each round builds a mutation->edges index, finds all candidate groups of k >= 2 siblings
/// whose mutation sets share a non-empty intersection, selects a disjoint set of groups
/// (no edge in two groups), and merges each group under one new internal node. The number
/// of children drops by at least one per round, so the loop always terminates.
///
/// Groups with k > 2 edges are handled directly in one merge, avoiding the chain of binary
/// internal nodes that would result from repeated pairwise merging.
///
/// Returns number of new internal nodes created.
fn merge_single_polytomy(
  graph: &mut GraphAncestral,
  partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
  node_key: GraphNodeKey,
) -> Result<usize, Report> {
  let mut nodes_created = 0;

  loop {
    let child_edges = collect_child_edge_keys(graph, node_key);
    if child_edges.len() <= 2 {
      break;
    }

    let index = build_mutation_index(partitions, &child_edges);
    let candidates = find_merge_groups(&index, partitions.len());
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
      merge_sibling_group(graph, partitions, node_key, &group)?;
      nodes_created += 1;
    }
  }

  Ok(nodes_created)
}

/// Collect outbound edge keys for a node.
fn collect_child_edge_keys(graph: &GraphAncestral, node_key: GraphNodeKey) -> Vec<GraphEdgeKey> {
  let node = graph.get_node(node_key).expect("Node must exist");
  let node = node.read_arc();
  node.outbound().to_vec()
}

/// Unified key for the mutation index: a substitution or an exact-match indel.
///
/// Substitutions are keyed by their full (ref, pos, qry) triple. Indels are keyed by their
/// exact (range, seq, deletion) value -- only identical indels count as shared. Partial
/// overlaps are ignored; this is conservative but cheap.
#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
enum MutationKey {
  Sub(Sub),
  InDel(InDel),
}

/// A group of k >= 2 sibling edges whose mutation sets share a non-empty intersection.
struct MergeGroup {
  /// Sorted, deduplicated edge keys.
  edges: Vec<GraphEdgeKey>,
  /// Shared substitutions per partition (index matches `partitions` slice).
  shared_subs: Vec<Vec<Sub>>,
  /// Shared indels per partition.
  shared_indels: Vec<Vec<InDel>>,
  /// Total shared mutations across all partitions.
  total_shared: usize,
}

/// Maps `(partition_index, MutationKey)` to the list of child edges carrying that mutation.
type MutationIndex = BTreeMap<(usize, MutationKey), Vec<GraphEdgeKey>>;

/// Build the mutation index for the given set of child edges.
///
/// Each substitution and each indel on each edge is entered under its (partition, key) bucket.
/// Buckets with a single entry contribute no pair scores; buckets with >= 2 entries contribute
/// +1 to all pairwise scores within the bucket.
fn build_mutation_index(
  partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
  child_edges: &[GraphEdgeKey],
) -> MutationIndex {
  let mut index: MutationIndex = BTreeMap::new();
  for &edge_key in child_edges {
    for (pi, partition) in partitions.iter().enumerate() {
      let partition = partition.read_arc();
      let empty_subs: &[Sub] = &[];
      let edge = partition.edges.get(&edge_key);
      for sub in edge.map_or(empty_subs, |e| e.fitch_subs()) {
        index.entry((pi, MutationKey::Sub(sub.clone()))).or_default().push(edge_key);
      }
      let empty_indels: &[InDel] = &[];
      for indel in edge.map_or(empty_indels, |e| e.indels.as_slice()) {
        index.entry((pi, MutationKey::InDel(indel.clone()))).or_default().push(edge_key);
      }
    }
  }
  index
}

/// Find all candidate merge groups from the mutation index.
///
/// For each index bucket with k >= 2 edges, the full set of edges in that bucket is a
/// candidate group. Duplicate candidate sets (same set of edges, different mutation) are
/// collapsed into one by using the sorted edge list as a dedup key. Each surviving candidate
/// is scored by intersecting all (partition, mutation) buckets that cover every edge in the
/// group.
fn find_merge_groups(index: &MutationIndex, n_partitions: usize) -> Vec<MergeGroup> {
  let mut seen: HashSet<Vec<GraphEdgeKey>> = HashSet::new();
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
    let (shared_subs, shared_indels) = shared_mutations_for_group(index, n_partitions, &group_edges);
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

/// Collect per-partition shared substitutions and indels for a group of edges.
///
/// A mutation is shared by the group if every edge in `group` appears in the bucket for
/// that (partition, mutation) key.
fn shared_mutations_for_group(
  index: &MutationIndex,
  n_partitions: usize,
  group: &[GraphEdgeKey],
) -> (Vec<Vec<Sub>>, Vec<Vec<InDel>>) {
  let mut shared_subs: Vec<Vec<Sub>> = vec![Vec::new(); n_partitions];
  let mut shared_indels: Vec<Vec<InDel>> = vec![Vec::new(); n_partitions];
  for ((pi, key), edges) in index {
    if group.iter().all(|e| edges.contains(e)) {
      match key {
        MutationKey::Sub(sub) => shared_subs[*pi].push(sub.clone()),
        MutationKey::InDel(indel) => shared_indels[*pi].push(indel.clone()),
      }
    }
  }
  (shared_subs, shared_indels)
}

/// Select a maximal disjoint set of merge groups in one round.
///
/// Sorts candidates by descending `total_shared`, breaking ties by descending group size.
/// Greedily picks groups where no edge has been claimed yet.
fn greedy_disjoint_group_matching(mut groups: Vec<MergeGroup>) -> Vec<MergeGroup> {
  groups.sort_unstable_by(|a, b| {
    b.total_shared
      .cmp(&a.total_shared)
      .then(b.edges.len().cmp(&a.edges.len()))
  });
  let mut used: HashSet<GraphEdgeKey> = HashSet::new();
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

/// Remaining mutations for one child edge after removing the shared mutations.
struct ChildEdgeData {
  remaining_subs: Vec<Sub>,
  remaining_indels: Vec<InDel>,
}

/// Merge k >= 2 siblings under a new internal node.
///
/// Creates a new node N between parent P and children C_0 ... C_{k-1}:
/// - Edge P->N carries the shared mutations.
/// - Edges N->C_i carry only the remaining unique mutations.
/// - Branch length of P->N is the Jukes-Cantor corrected distance from
///   `total_shared / total_alignment_length`.
/// - Branch lengths of N->C_i = JC(remaining_i / total_alignment_length).
///
/// # Model assumption
///
/// The `prune` command initialises all partitions with JC69 (`run_prune()`),
/// so applying the Jukes-Cantor correction is exact. When this function is
/// reused from `optimize` under a non-JC69 model, JC69 remains a strictly
/// better approximation than the raw p-distance because any symmetric
/// substitution process underestimates true evolutionary distance by the
/// same mechanism (back-mutations and parallel substitutions).
fn merge_sibling_group(
  graph: &mut GraphAncestral,
  partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
  parent_key: GraphNodeKey,
  group: &MergeGroup,
) -> Result<(), Report> {
  let child_keys: Vec<GraphNodeKey> = group
    .edges
    .iter()
    .map(|&ek| graph.get_target_node_key(ek))
    .collect::<Result<Vec<_>, _>>()?;

  let per_edge_data: Vec<Vec<ChildEdgeData>> = partitions
    .iter()
    .enumerate()
    .map(|(pi, partition)| {
      let partition = partition.read_arc();
      let empty_subs: &[Sub] = &[];
      let empty_indels: &[InDel] = &[];
      group
        .edges
        .iter()
        .map(|&ek| {
          let edge = partition.edges.get(&ek);
          let all_subs = edge.map_or(empty_subs, |e| e.fitch_subs());
          let all_indels = edge.map_or(empty_indels, |e| e.indels.as_slice());
          let shared_subs = &group.shared_subs[pi];
          let shared_indels = &group.shared_indels[pi];
          let remaining_subs = iterator_difference(all_subs, shared_subs).cloned().collect_vec();
          let remaining_indels = iterator_difference(all_indels, shared_indels.as_slice()).cloned().collect_vec();
          ChildEdgeData {
            remaining_subs,
            remaining_indels,
          }
        })
        .collect()
    })
    .collect();

  let total_alignment_length: usize = partitions.iter().map(|p| p.read_arc().length).sum();
  let n_states = partitions[0].read_arc().alphabet.n_canonical();
  let jc_bl = |count: usize| -> f64 {
    if total_alignment_length > 0 {
      jukes_cantor_distance(count as f64 / total_alignment_length as f64, n_states)
    } else {
      0.0
    }
  };

  let new_edge_bl = jc_bl(group.total_shared);
  let child_bls: Vec<f64> = (0..group.edges.len())
    .map(|i| {
      let remaining: usize = per_edge_data
        .iter()
        .map(|pi_data| pi_data[i].remaining_subs.len() + pi_data[i].remaining_indels.len())
        .sum();
      jc_bl(remaining)
    })
    .collect();

  for &ek in &group.edges {
    graph.remove_edge(ek)?;
  }

  let new_node_key = graph.add_node(NodeAncestral::default());

  let new_parent_edge_key = graph.add_edge(
    parent_key,
    new_node_key,
    EdgeAncestral {
      branch_length: Some(new_edge_bl),
    },
  )?;

  let new_child_edge_keys: Vec<GraphEdgeKey> = child_keys
    .iter()
    .zip(child_bls.iter())
    .map(|(&ck, &bl)| {
      graph.add_edge(
        new_node_key,
        ck,
        EdgeAncestral {
          branch_length: Some(bl),
        },
      )
    })
    .collect::<Result<Vec<_>, _>>()?;

  for (pi, partition_arc) in partitions.iter().enumerate() {
    let mut partition = partition_arc.write_arc();

    let mut new_node = SparseNodePartition::empty(&partition.alphabet);
    let parent_comp = partition.nodes[&parent_key].seq.composition.clone();
    new_node.seq.composition = parent_comp.clone();
    new_node.seq.fitch.composition = parent_comp;
    partition.nodes.entry(new_node_key).or_insert(new_node);

    for &ek in &group.edges {
      partition.edges.remove(&ek);
    }

    let parent_edge = partition.edges.entry(new_parent_edge_key).or_default();
    parent_edge.set_fitch_subs(group.shared_subs[pi].clone());
    parent_edge.indels = group.shared_indels[pi].clone();

    for (i, &new_ek) in new_child_edge_keys.iter().enumerate() {
      let child_edge = partition.edges.entry(new_ek).or_default();
      child_edge.set_fitch_subs(per_edge_data[pi][i].remaining_subs.clone());
      child_edge.indels = per_edge_data[pi][i].remaining_indels.clone();
    }
  }

  Ok(())
}
