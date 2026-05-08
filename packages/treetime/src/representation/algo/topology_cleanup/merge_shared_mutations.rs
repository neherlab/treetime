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
use std::sync::Arc;
use treetime_graph::edge::{GraphEdgeKey, HasBranchLength};
use treetime_graph::node::GraphNodeKey;
use treetime_utils::iterator::difference::iterator_difference;
use treetime_utils::iterator::intersection::iterator_intersection;

/// Merge sibling branches in polytomies that share identical substitutions.
///
/// Tree builders produce arbitrary binary resolutions of polytomies. When two sibling
/// branches carry identical substitutions, they represent redundant resolutions. This
/// function groups such siblings under a new internal node:
///
/// Before: P-A (subs: {G100T, A200C}), P-B (subs: {G100T, A200C, T300G}), P-C
/// After:  P-N (subs: {G100T, A200C}), N-A (subs: {}), N-B (subs: {T300G}), P-C
///
/// Greedily merges the pair with the most shared mutations, repeating until no pair
/// shares any mutations. Returns the total number of new internal nodes created.
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
    debug!("Merged {total_merged} sibling pairs sharing mutations");
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

/// Resolve one polytomy by greedy pairwise merging of siblings sharing mutations.
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

    let best = find_best_shared_mutation_pair(graph, partitions, &child_edges)?;
    let Some(best) = best else { break };

    debug!(
      "Merging siblings under node {node_key}: edges {} and {} share {} mutations",
      best.edge_key_a, best.edge_key_b, best.total_shared
    );

    merge_sibling_pair(graph, partitions, node_key, &best)?;
    nodes_created += 1;
  }

  Ok(nodes_created)
}

/// Collect outbound edge keys for a node.
fn collect_child_edge_keys(graph: &GraphAncestral, node_key: GraphNodeKey) -> Vec<GraphEdgeKey> {
  let node = graph.get_node(node_key).expect("Node must exist");
  let node = node.read_arc();
  node.outbound().to_vec()
}

/// Result of finding the best pair of siblings to merge.
struct SharedMutationPair {
  edge_key_a: GraphEdgeKey,
  edge_key_b: GraphEdgeKey,
  /// Shared substitutions per partition (indexed by partition position).
  shared_subs: Vec<Vec<Sub>>,
  total_shared: usize,
}

/// Find the pair of child edges with the most shared mutations across all partitions.
/// Returns None if no pair shares any mutations.
fn find_best_shared_mutation_pair(
  graph: &GraphAncestral,
  partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
  child_edges: &[GraphEdgeKey],
) -> Result<Option<SharedMutationPair>, Report> {
  let n = child_edges.len();
  let mut best: Option<SharedMutationPair> = None;

  for i in 0..n {
    for j in (i + 1)..n {
      let shared_subs = compute_shared_subs_across_partitions(partitions, child_edges[i], child_edges[j])?;
      let total_shared: usize = shared_subs.iter().map(Vec::len).sum();

      if total_shared > 0 {
        let is_better = best.as_ref().is_none_or(|b| total_shared > b.total_shared);
        if is_better {
          best = Some(SharedMutationPair {
            edge_key_a: child_edges[i],
            edge_key_b: child_edges[j],
            shared_subs,
            total_shared,
          });
        }
      }
    }
  }

  Ok(best)
}

/// Compute the intersection of substitutions on two edges, per partition.
fn compute_shared_subs_across_partitions(
  partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
  edge_key_a: GraphEdgeKey,
  edge_key_b: GraphEdgeKey,
) -> Result<Vec<Vec<Sub>>, Report> {
  let empty: &[Sub] = &[];
  partitions
    .iter()
    .map(|partition| {
      let partition = partition.read_arc();
      let subs_a = partition.edges.get(&edge_key_a).map_or(empty, |e| e.fitch_subs());
      let subs_b = partition.edges.get(&edge_key_b).map_or(empty, |e| e.fitch_subs());
      Ok(iterator_intersection(subs_a, subs_b).cloned().collect_vec())
    })
    .collect()
}

/// Saved partition data for two child edges during merge.
struct ChildPartitionData {
  remaining_subs_a: Vec<Sub>,
  remaining_subs_b: Vec<Sub>,
  indels_a: Vec<InDel>,
  indels_b: Vec<InDel>,
}

/// Merge two siblings under a new internal node.
///
/// Creates a new node N between parent P and children A, B:
/// - Edge P-N carries the shared mutations
/// - Edges N-A and N-B carry only the remaining unique mutations
/// - Branch length of P-N is the Jukes-Cantor corrected evolutionary
///   distance from the pooled p-distance `total_shared_mutations / total_alignment_length`
///   (see [`jukes_cantor_distance`])
/// - Branch lengths of N-A, N-B = max(0, original - new_edge_bl)
///
/// # Model assumption
///
/// The `prune` command initialises all partitions with JC69 (`run_prune()`),
/// so applying the Jukes-Cantor correction is exact. When this function is
/// reused from `optimize` under a non-JC69 model, JC69 remains a strictly
/// better approximation than the raw p-distance because any symmetric
/// substitution process underestimates true evolutionary distance by the
/// same mechanism (back-mutations and parallel substitutions).
fn merge_sibling_pair(
  graph: &mut GraphAncestral,
  partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
  parent_key: GraphNodeKey,
  pair: &SharedMutationPair,
) -> Result<(), Report> {
  // Read actual child keys from graph edges
  let child_key_a = graph.get_target_node_key(pair.edge_key_a)?;
  let child_key_b = graph.get_target_node_key(pair.edge_key_b)?;

  // Read original branch lengths before removing edges
  let bl_a = graph
    .get_edge(pair.edge_key_a)
    .and_then(|e| e.read_arc().payload().read_arc().branch_length());
  let bl_b = graph
    .get_edge(pair.edge_key_b)
    .and_then(|e| e.read_arc().payload().read_arc().branch_length());

  let old_partition_data: Vec<_> = partitions
    .iter()
    .enumerate()
    .map(|(pi, partition)| {
      let partition = partition.read_arc();
      let empty: &[Sub] = &[];
      let edge_a = partition.edges.get(&pair.edge_key_a);
      let edge_b = partition.edges.get(&pair.edge_key_b);
      let subs_a = edge_a.map_or(empty, |e| e.fitch_subs());
      let subs_b = edge_b.map_or(empty, |e| e.fitch_subs());
      let indels_a = edge_a.map(|e| e.indels.clone()).unwrap_or_default();
      let indels_b = edge_b.map(|e| e.indels.clone()).unwrap_or_default();
      let shared = &pair.shared_subs[pi];
      let remaining_a = iterator_difference(subs_a, shared).cloned().collect_vec();
      let remaining_b = iterator_difference(subs_b, shared).cloned().collect_vec();
      Ok(ChildPartitionData {
        remaining_subs_a: remaining_a,
        remaining_subs_b: remaining_b,
        indels_a,
        indels_b,
      })
    })
    .collect::<Result<Vec<_>, Report>>()?;

  // Compute new branch length from the pooled p-distance across partitions.
  // The raw Hamming ratio shared/length underestimates evolutionary distance
  // because multiple substitutions at the same site can revert or mask earlier
  // ones. The Jukes-Cantor 1969 correction inverts the JC substitution model
  // to recover expected substitutions per site from observed differences.
  //
  // Pooling across partitions assumes they share an alphabet size, which holds
  // for every current caller: `prune` hardcodes JC69 over a single alphabet,
  // and `optimize` reuses this path with a single partition set under one
  // model name. `n_states` is read from the first partition accordingly.
  let total_alignment_length: usize = partitions.iter().map(|p| p.read_arc().length).sum();
  let new_edge_bl = if total_alignment_length > 0 {
    let p = pair.total_shared as f64 / total_alignment_length as f64;
    let n_states = partitions[0].read_arc().alphabet.n_canonical();
    jukes_cantor_distance(p, n_states)
  } else {
    0.0
  };

  // Adjust children branch lengths
  let bl_a_adjusted = bl_a.map(|bl| f64::max(0.0, bl - new_edge_bl));
  let bl_b_adjusted = bl_b.map(|bl| f64::max(0.0, bl - new_edge_bl));

  // Remove old edges from parent to children
  graph.remove_edge(pair.edge_key_a)?;
  graph.remove_edge(pair.edge_key_b)?;

  // Create new internal node
  let new_node_key = graph.add_node(NodeAncestral::default());

  // Edge parent-new_node carries shared mutations
  let new_parent_edge_key = graph.add_edge(
    parent_key,
    new_node_key,
    EdgeAncestral {
      branch_length: Some(new_edge_bl),
    },
  )?;

  // Edge new_node-child_a carries remaining mutations
  let new_edge_a_key = graph.add_edge(
    new_node_key,
    child_key_a,
    EdgeAncestral {
      branch_length: bl_a_adjusted,
    },
  )?;

  // Edge new_node-child_b carries remaining mutations
  let new_edge_b_key = graph.add_edge(
    new_node_key,
    child_key_b,
    EdgeAncestral {
      branch_length: bl_b_adjusted,
    },
  )?;

  // Update partition node and edge data
  for (pi, partition) in partitions.iter().enumerate() {
    let mut partition = partition.write_arc();
    let data = &old_partition_data[pi];

    // Add node entry for the new internal node. Copy the parent's composition
    // so that the backward pass in combine_messages() computes correct fixed-site
    // contributions. An empty composition would zero out ~99% of the log-likelihood
    // signal for the merge-created node and all ancestors.
    let mut new_node = SparseNodePartition::empty(&partition.alphabet);
    let parent_comp = partition.nodes[&parent_key].seq.composition.clone();
    new_node.seq.composition = parent_comp.clone();
    new_node.seq.fitch.composition = parent_comp;
    partition.nodes.entry(new_node_key).or_insert(new_node);

    // Remove old edge entries
    partition.edges.remove(&pair.edge_key_a);
    partition.edges.remove(&pair.edge_key_b);

    // Insert shared mutations on parent-to-new-node edge (no indels - only point mutations are shared)
    let parent_edge = partition.edges.entry(new_parent_edge_key).or_default();
    parent_edge.set_fitch_subs(pair.shared_subs[pi].clone());

    // Insert remaining mutations and preserved indels on new-node-to-child edges
    let edge_a = partition.edges.entry(new_edge_a_key).or_default();
    edge_a.set_fitch_subs(data.remaining_subs_a.clone());
    edge_a.indels = data.indels_a.clone();

    let edge_b = partition.edges.entry(new_edge_b_key).or_default();
    edge_b.set_fitch_subs(data.remaining_subs_b.clone());
    edge_b.indels = data.indels_b.clone();
  }

  Ok(())
}
