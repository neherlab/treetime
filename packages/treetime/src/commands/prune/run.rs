use crate::alphabet::alphabet::Alphabet;
use crate::commands::prune::args::TreetimePruneArgs;
use crate::gtr::get_gtr::{GtrModelName, log_gtr, write_gtr_json};
use crate::make_error;
use crate::representation::algo::topology_cleanup::collapse::collapse_edge;
use crate::representation::algo::topology_cleanup::merge_shared_mutations::merge_shared_mutation_branches;
use crate::representation::partition::fitch::PartitionFitch;
use crate::representation::partition::marginal_dense::PartitionMarginalDense;
use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
use crate::representation::payload::ancestral::GraphAncestral;
use eyre::Report;
use itertools::Itertools;
use log::debug;
use maplit::btreeset;
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

    let fitch = PartitionFitch::compress(&graph, 0, alphabet, &aln)?;
    let gtr = fitch.resolve_gtr(&graph, GtrModelName::JC69)?;
    log_gtr(&gtr, GtrModelName::JC69);
    write_gtr_json(&gtr, GtrModelName::JC69, outdir, None)?;
    let partition = fitch.into_marginal_sparse(gtr, &graph)?;
    vec![Arc::new(RwLock::new(partition))]
  } else {
    vec![]
  };

  let node_names = parse_node_names(
    prune_nodes_list.as_ref(),
    *prune_nodes_list_delimiter,
    prune_nodes_list_file.as_ref(),
    *prune_nodes_list_file_delimiter,
  )?;

  // Prune first: collapsing short/empty branches creates larger polytomies,
  // exposing more siblings for shared-mutation merging.
  prune_nodes(&mut graph, &partitions, *prune_short, *prune_empty, &node_names)?;

  if *merge_shared_mutations {
    merge_shared_mutation_branches(&mut graph, &partitions)?;
    graph.build()?;
  }

  write_graph(outdir, &graph)?;

  Ok(())
}

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
pub(super) fn get_edge_num_muts(
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

pub(super) fn collapse_sparse_edges_from_leaf_recursive(
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

fn should_collapse_parent(graph: &GraphAncestral, node_key: GraphNodeKey) -> bool {
  graph.has_at_most_one_child(node_key) && !graph.is_root(node_key)
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
