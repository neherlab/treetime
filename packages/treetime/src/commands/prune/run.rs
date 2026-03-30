use crate::alphabet::alphabet::Alphabet;
use crate::commands::ancestral::fitch::{compress_sequences, get_common_length};
use crate::commands::prune::args::TreetimePruneArgs;
use crate::gtr::get_gtr::{GtrModelName, JC69Params, get_gtr_sparse, jc69, write_gtr_json};
use crate::make_error;
use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
use crate::representation::partition::traits::PartitionBranchOps;
use crate::representation::payload::ancestral::{EdgeAncestral, GraphAncestral, NodeAncestral};
use crate::representation::payload::sparse::SparseNodePartition;
use crate::seq::indel::InDel;
use crate::seq::mutation::{Sub, compose_substitutions};
use eyre::Report;
use itertools::Itertools;
use log::debug;
use maplit::{btreemap, btreeset};
use parking_lot::RwLock;
use serde::Serialize;
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use treetime_graph::edge::{GraphEdge, GraphEdgeKey, HasBranchLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, GraphNodeKey, Named};
use treetime_io::fasta::read_many_fasta;
use treetime_io::nex::{NexWriteOptions, nex_write_file};
use treetime_io::nwk::{EdgeToNwk, NodeToNwk, NwkWriteOptions, nwk_read_file, nwk_write_file};
use treetime_io::parse_delimited::{parse_delimited_file, parse_delimited_str};
use treetime_utils::iterator::difference::iterator_difference;
use treetime_utils::iterator::intersection::iterator_intersection;

fn validate_args(args: &TreetimePruneArgs) -> Result<(), Report> {
  if args.prune_empty && args.input_fastas.is_empty() {
    return make_error!(
      "The --prune-empty requires --aln. Without sequence data, it's not possible to determine which branches lack mutations."
    );
  }

  if args.merge_shared_mutations && args.input_fastas.is_empty() {
    return make_error!(
      "The --merge-shared-mutations requires --aln. Without sequence data, it's not possible to determine which branches share mutations."
    );
  }

  Ok(())
}

pub fn run_prune(args: &TreetimePruneArgs) -> Result<(), Report> {
  validate_args(args)?;

  let TreetimePruneArgs {
    input_fastas,
    tree,
    alphabet,
    outdir,
    prune_short,
    prune_empty,
    merge_shared_mutations,
    prune_nodes_list,
    prune_nodes_list_delimiter,
    prune_nodes_list_file,
    prune_nodes_list_file_delimiter,
  } = args;

  let mut graph: GraphAncestral = nwk_read_file(tree)?;

  let needs_sequences = *prune_empty || *merge_shared_mutations;
  let partitions: Vec<Arc<RwLock<PartitionMarginalSparse>>> = if needs_sequences {
    let alphabet = Alphabet::new(alphabet.unwrap_or_default())?;
    let aln = read_many_fasta(input_fastas, &alphabet)?;

    let partitions = vec![Arc::new(RwLock::new(PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?, // FIXME: dummy temporary gtr should not be needed here
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    compress_sequences(&graph, &partitions, &aln)?;

    // FIXME: chicken & egg problem: to get a gtr we need partitions, to get partitions we need a gtr
    // FIXME: spaghetti code: dummy gtr is replaced by real gtr here
    for partition in &partitions {
      let gtr = get_gtr_sparse(&GtrModelName::JC69, partition, &graph)?;
      write_gtr_json(&gtr, GtrModelName::JC69, outdir, None)?;
      partition.write_arc().gtr = gtr;
    }

    partitions
  } else {
    vec![]
  };

  let node_names = parse_node_names(
    prune_nodes_list.as_ref(),
    *prune_nodes_list_delimiter,
    prune_nodes_list_file.as_ref(),
    *prune_nodes_list_file_delimiter,
  )?;

  if *merge_shared_mutations {
    merge_shared_mutation_branches(&mut graph, &partitions)?;
    graph.build()?;
  }

  prune_nodes(&mut graph, &partitions, *prune_short, *prune_empty, &node_names)?;

  write_graph(outdir, &graph)?;

  Ok(())
}

fn parse_node_names(
  prune_nodes_list: Option<&String>,
  prune_nodes_list_delimiter: char,
  prune_nodes_list_file: Option<&PathBuf>,
  prune_nodes_list_file_delimiter: char,
) -> Result<BTreeSet<String>, Report> {
  let mut node_names = btreeset! {};

  if let Some(prune_nodes_list) = prune_nodes_list {
    let names: Vec<String> = parse_delimited_str(prune_nodes_list, prune_nodes_list_delimiter as u8).try_collect()?;
    node_names.extend(names);
  }

  if let Some(prune_nodes_list_file) = prune_nodes_list_file {
    let names: Vec<String> =
      parse_delimited_file(prune_nodes_list_file, prune_nodes_list_file_delimiter as u8)?.try_collect()?;
    node_names.extend(names);
  }

  Ok(node_names)
}

pub(super) fn prune_nodes(
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
///
/// Uses `PartitionBranchOps::edge_subs()` to derive mutations from the
/// current partition state rather than reading stale Fitch-era `subs`.
pub(super) fn get_edge_num_muts(
  graph: &GraphAncestral,
  partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
  edge_key: GraphEdgeKey,
) -> Result<Option<usize>, Report> {
  let mut total_muts = 0;
  let mut found_any = false;

  for partition in partitions {
    let partition = partition.read_arc();
    if partition.edges.contains_key(&edge_key) {
      total_muts += partition.edge_subs(graph, edge_key)?.len();
      found_any = true;
    }
  }

  Ok(found_any.then_some(total_muts))
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

      // Skip leaves in this pass
      if target_is_leaf {
        return Ok(None);
      }

      let weight = edge.payload().read_arc().branch_length();
      let should_prune_short = matches!((prune_short, weight), (Some(threshold), Some(weight)) if weight < threshold);

      let should_prune_empty = prune_empty && get_edge_num_muts(graph, partitions, edge.key())? == Some(0);

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

  edges_to_collapse.into_iter().try_for_each(|edge_key| {
    debug!("Collapsing internal edge: {edge_key}");
    collapse_sparse_edge(graph, partitions, edge_key)
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

pub(super) fn collapse_sparse_edges_from_leaf_recursive(
  graph: &mut GraphAncestral,
  partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
  edge_key: GraphEdgeKey,
) -> Result<(), Report> {
  let mut current_edge_key = edge_key;

  loop {
    let parent_node_key = graph.get_source_node_key(current_edge_key)?;

    collapse_sparse_edge(graph, partitions, current_edge_key)?;

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

fn should_collapse_parent(graph: &GraphAncestral, node_key: GraphNodeKey) -> bool {
  graph.has_at_most_one_child(node_key) && !graph.is_root(node_key)
}

fn collapse_sparse_edge(
  graph: &mut GraphAncestral,
  partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
  edge_key: GraphEdgeKey,
) -> Result<(), Report> {
  let (_, removed_edge, new_edges) = graph.collapse_edge(edge_key)?;
  let removed_edge = removed_edge.payload().read_arc();

  for new_edge in new_edges {
    let new_edge_key = new_edge.read_arc().key();
    let mut new_edge = new_edge.write_arc().payload().write_arc();

    // Sum branch lengths
    if let (Some(bl1), Some(bl2)) = (removed_edge.branch_length(), new_edge.branch_length()) {
      new_edge.set_branch_length(Some(bl1 + bl2));
    }

    // Compose substitutions and merge indels: parent edge (removed) + child edge
    for partition in partitions {
      let mut partition = partition.write_arc();
      let removed_edge_data = partition.edges[&edge_key].clone();
      let child_edge = partition.edges.entry(new_edge_key).or_default();
      child_edge.subs = compose_substitutions(&removed_edge_data.subs, &child_edge.subs)?;
      let mut merged_indels = removed_edge_data.indels;
      merged_indels.append(&mut child_edge.indels);
      child_edge.indels = merged_indels;
    }
  }

  Ok(())
}

/// Merge sibling branches in polytomies that share identical substitutions.
///
/// Tree builders produce arbitrary binary resolutions of polytomies. When two sibling
/// branches carry identical substitutions, they represent redundant resolutions. This
/// function groups such siblings under a new internal node:
///
/// Before: P → A (subs: {G100T, A200C}), P → B (subs: {G100T, A200C, T300G}), P → C
/// After:  P → N (subs: {G100T, A200C}), N → A (subs: {}), N → B (subs: {T300G}), P → C
///
/// Greedily merges the pair with the most shared mutations, repeating until no pair
/// shares any mutations. Returns the total number of new internal nodes created.
pub(super) fn merge_shared_mutation_branches(
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
      let shared_subs = compute_shared_subs_across_partitions(graph, partitions, child_edges[i], child_edges[j])?;
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
///
/// Uses `PartitionBranchOps::edge_subs()` to derive mutations from the
/// current partition state rather than reading stale Fitch-era `subs`.
fn compute_shared_subs_across_partitions(
  graph: &GraphAncestral,
  partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
  edge_key_a: GraphEdgeKey,
  edge_key_b: GraphEdgeKey,
) -> Result<Vec<Vec<Sub>>, Report> {
  partitions
    .iter()
    .map(|partition| {
      let partition = partition.read_arc();
      let subs_a = partition.edge_subs(graph, edge_key_a)?;
      let subs_b = partition.edge_subs(graph, edge_key_b)?;
      Ok(iterator_intersection(&subs_a, &subs_b).cloned().collect_vec())
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
/// - Edge P → N carries the shared mutations
/// - Edges N → A and N → B carry only the remaining unique mutations
/// - Branch length of P → N = total_shared_mutations / total_alignment_length
/// - Branch lengths of N → A, N → B = max(0, original - new_edge_bl)
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

  // Save original partition data for both edges before removal.
  // Uses edge_subs() to derive mutations from current state rather than
  // reading stale Fitch-era subs directly.
  let old_partition_data: Vec<_> = partitions
    .iter()
    .enumerate()
    .map(|(pi, partition)| {
      let partition = partition.read_arc();
      let edge_a = partition.edges.get(&pair.edge_key_a);
      let edge_b = partition.edges.get(&pair.edge_key_b);
      let subs_a = partition.edge_subs(graph, pair.edge_key_a)?;
      let subs_b = partition.edge_subs(graph, pair.edge_key_b)?;
      let indels_a = edge_a.map(|e| e.indels.clone()).unwrap_or_default();
      let indels_b = edge_b.map(|e| e.indels.clone()).unwrap_or_default();
      let shared = &pair.shared_subs[pi];
      let remaining_a = iterator_difference(&subs_a, shared).cloned().collect_vec();
      let remaining_b = iterator_difference(&subs_b, shared).cloned().collect_vec();
      Ok(ChildPartitionData {
        remaining_subs_a: remaining_a,
        remaining_subs_b: remaining_b,
        indels_a,
        indels_b,
      })
    })
    .collect::<Result<Vec<_>, Report>>()?;

  // Compute new branch length: shared_mutations / alignment_length
  let total_alignment_length: usize = partitions.iter().map(|p| p.read_arc().length).sum();
  let new_edge_bl = if total_alignment_length > 0 {
    pair.total_shared as f64 / total_alignment_length as f64
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

  // Create edge: parent → new_node (carries shared mutations)
  let new_parent_edge_key = graph.add_edge(
    parent_key,
    new_node_key,
    EdgeAncestral {
      branch_length: Some(new_edge_bl),
    },
  )?;

  // Create edge: new_node → child_a (carries remaining mutations)
  let new_edge_a_key = graph.add_edge(
    new_node_key,
    child_key_a,
    EdgeAncestral {
      branch_length: bl_a_adjusted,
    },
  )?;

  // Create edge: new_node → child_b (carries remaining mutations)
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

    // Add node entry for the new internal node so edge_subs_from_graph()
    // can reconstruct states on edges connected to it.
    let empty_node = SparseNodePartition::empty(&partition.alphabet);
    partition.nodes.entry(new_node_key).or_insert(empty_node);

    // Remove old edge entries
    partition.edges.remove(&pair.edge_key_a);
    partition.edges.remove(&pair.edge_key_b);

    // Insert shared mutations on parent-to-new-node edge (no indels - only point mutations are shared)
    let parent_edge = partition.edges.entry(new_parent_edge_key).or_default();
    parent_edge.subs = pair.shared_subs[pi].clone();

    // Insert remaining mutations and preserved indels on new-node-to-child edges
    let edge_a = partition.edges.entry(new_edge_a_key).or_default();
    edge_a.subs = data.remaining_subs_a.clone();
    edge_a.indels = data.indels_a.clone();

    let edge_b = partition.edges.entry(new_edge_b_key).or_default();
    edge_b.subs = data.remaining_subs_b.clone();
    edge_b.indels = data.indels_b.clone();
  }

  Ok(())
}

fn write_graph<N, E, D>(outdir: impl AsRef<Path>, graph: &Graph<N, E, D>) -> Result<(), Report>
where
  N: GraphNode + NodeToNwk + Serialize,
  E: GraphEdge + EdgeToNwk + Serialize,
  D: Send + Sync + Default + Serialize,
{
  nwk_write_file(
    outdir.as_ref().join("pruned_tree.nwk"),
    graph,
    &NwkWriteOptions::default(),
  )?;

  nex_write_file(
    outdir.as_ref().join("pruned_tree.nexus"),
    graph,
    &NexWriteOptions::default(),
  )?;

  Ok(())
}
