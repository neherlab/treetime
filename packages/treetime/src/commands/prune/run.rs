use crate::alphabet::alphabet::Alphabet;
use crate::commands::ancestral::fitch::{compress_sequences, get_common_length};
use crate::commands::prune::args::TreetimePruneArgs;
use crate::graph::edge::{GraphEdge, GraphEdgeKey, HasBranchLength};
use crate::graph::graph::Graph;
use crate::graph::node::{GraphNode, GraphNodeKey, Named};
use crate::gtr::get_gtr::{GtrModelName, JC69Params, get_gtr, jc69};
use crate::io::fasta::read_many_fasta;
use crate::io::nex::{NexWriteOptions, nex_write_file};
use crate::io::nwk::{EdgeToNwk, NodeToNwk, NwkWriteOptions, nwk_read_file, nwk_write_file};
use crate::io::parse_delimited::{parse_delimited_file, parse_delimited_str};
use crate::make_error;
use crate::representation::graph_ancestral::GraphAncestral;
use crate::representation::partition_marginal_sparse::PartitionMarginalSparse;
use eyre::Report;
use itertools::Itertools;
use log::debug;
use maplit::{btreemap, btreeset};
use parking_lot::RwLock;
use serde::Serialize;
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use treetime_utils::iterator::union::iterator_union;

fn validate_args(args: &TreetimePruneArgs) -> Result<(), Report> {
  if args.prune_empty && args.input_fastas.is_empty() {
    return make_error!(
      "The --prune-empty requires --aln. Without sequence data, it's not possible to determine which branches lack mutations."
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
    prune_nodes_list,
    prune_nodes_list_delimiter,
    prune_nodes_list_file,
    prune_nodes_list_file_delimiter,
  } = args;

  let mut graph: GraphAncestral = nwk_read_file(tree)?;

  let partitions: Vec<Arc<RwLock<PartitionMarginalSparse>>> = if *prune_empty {
    assert!(
      !input_fastas.is_empty(),
      "Input FASTA files are required for pruning empty branches."
    );
    let alphabet = Alphabet::new(alphabet.unwrap_or_default(), false)?;
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
      let gtr = get_gtr(&GtrModelName::JC69, partition, &graph)?;
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
    node_names.extend(parse_delimited_str(prune_nodes_list, prune_nodes_list_delimiter as u8));
  }

  if let Some(prune_nodes_list_file) = prune_nodes_list_file {
    node_names.extend(parse_delimited_file(
      prune_nodes_list_file,
      prune_nodes_list_file_delimiter as u8,
    )?);
  }

  Ok(node_names)
}

pub(super) fn prune_nodes(
  graph: &mut GraphAncestral,
  partitions: &Vec<Arc<RwLock<PartitionMarginalSparse>>>,
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

pub(super) fn get_edge_num_muts(
  partitions: &Vec<Arc<RwLock<PartitionMarginalSparse>>>,
  edge_key: GraphEdgeKey,
) -> Option<usize> {
  let mut total_muts = 0;
  let mut found_any = false;

  for partition in partitions {
    let partition = partition.read_arc();
    if let Some(edge_partition) = partition.edges.get(&edge_key) {
      total_muts += edge_partition.subs.len();
      found_any = true;
    }
  }

  found_any.then_some(total_muts)
}

fn prune_internal_nodes(
  graph: &mut GraphAncestral,
  partitions: &Vec<Arc<RwLock<PartitionMarginalSparse>>>,
  prune_short: Option<f64>,
  prune_empty: bool,
  node_names: &BTreeSet<String>,
) -> Result<(), Report> {
  #[allow(clippy::needless_collect)]
  let edges_to_collapse: Vec<_> = graph
    .get_edges()
    .iter()
    .filter_map(|edge| {
      let edge = edge.read_arc();
      let target_is_leaf = graph.is_leaf(edge.target());

      // Skip leaves in this pass
      if target_is_leaf {
        return None;
      }

      let weight = edge.payload().read_arc().branch_length();
      let should_prune_short = matches!((prune_short, weight), (Some(threshold), Some(weight)) if weight < threshold);

      let should_prune_empty = prune_empty && get_edge_num_muts(partitions, edge.key()) == Some(0);

      let target_node = graph.get_node(edge.target())?.read_arc().payload().read_arc();
      let name = target_node.name();
      let should_prune_by_name = name.is_some_and(|name| node_names.contains(name.as_ref()));

      let should_prune = should_prune_short || should_prune_empty || should_prune_by_name;
      should_prune.then(|| edge.key())
    })
    .collect();

  edges_to_collapse.into_iter().try_for_each(|edge_key| {
    debug!("Collapsing internal edge: {edge_key}");
    collapse_sparse_edge(graph, partitions, edge_key)
  })
}

fn prune_leaves(
  graph: &mut GraphAncestral,
  partitions: &Vec<Arc<RwLock<PartitionMarginalSparse>>>,
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
  partitions: &Vec<Arc<RwLock<PartitionMarginalSparse>>>,
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
  partitions: &Vec<Arc<RwLock<PartitionMarginalSparse>>>,
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

    // Union of substitutions per partition
    for partition in partitions {
      let mut partition = partition.write_arc();
      let removed_edge_subs = partition.edges[&edge_key].subs.clone(); // FIXME: avoid clone
      let new_edge = partition.edges.entry(new_edge_key).or_default();
      new_edge.subs = iterator_union(&removed_edge_subs, &new_edge.subs)
        .cloned()
        .collect_vec();
    }
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
