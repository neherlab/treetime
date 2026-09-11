use crate::alphabet::alphabet::Alphabet;
use crate::commands::prune::args::TreetimePruneArgs;
use crate::commands::prune::result::{EdgeOut, PruneGraphData, PruneNodeOut, PruneOutputMaps, PruneResult};
use crate::commands::shared::output::OutputSelection;
use crate::commands::shared::resolve_outputs::ResolveOutputs;
use crate::commands::shared::tree_output::write_prune_tree_outputs;
use crate::gtr::get_gtr::{GtrModelName, GtrOutput, write_gtr_json};
use crate::make_error;
use crate::partition::traits::PartitionBranchOps;
use crate::prune::pipeline::{self, PruneInput, PruneParams};
use crate::seq::mutation::MutationTrack;
use eyre::{Report, WrapErr};
use itertools::Itertools;
use log::warn;
use maplit::btreeset;
use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::path::PathBuf;
use std::sync::Arc;
use treetime_graph::edge::{GraphEdge, GraphEdgeKey};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, GraphNodeKey};
use treetime_io::fasta::read_many_fasta;
use treetime_io::graph::TreeWriteKind;
use treetime_io::nwk::CommentProviders;
use treetime_io::nwk::nwk_read_file;
use treetime_io::parse_delimited::{parse_delimited_file, parse_delimited_str};

use crate::payload::ancestral::GraphAncestral;

pub fn run_prune(
  args: &TreetimePruneArgs,
  progress: &dyn crate::progress::ProgressSink,
) -> Result<PruneResult, Report> {
  validate_args(args)?;

  progress.check_cancelled()?;
  progress.report("Reading input", 0.0, "");

  let parse = nwk_read_file(args.tree())?;
  let graph: GraphAncestral = parse.graph;
  let confidences = parse.confidences;
  let names = parse.names;
  let branch_lengths_input = parse.branch_lengths;
  let input_order = leaf_order(&graph, &names)?;
  let alphabet = Alphabet::new(args.alphabet_args.alphabet.unwrap_or_default())?;

  let resolved = args.resolve_outputs()?;

  let needs_sequences = args.prune_empty || args.merge_shared_mutations;
  let sequences = if needs_sequences && !args.alignment.alignment.is_empty() {
    Some(read_many_fasta(&args.alignment.alignment, &alphabet)?)
  } else {
    None
  };

  let node_names = parse_node_names(
    args.prune_nodes_list.as_ref(),
    args.prune_nodes_list_delimiter,
    args.prune_nodes_list_file.as_ref(),
    args.prune_nodes_list_file_delimiter,
  )?;

  let params = PruneParams {
    prune_short: args.prune_short,
    prune_empty: args.prune_empty,
    merge_shared_mutations: args.merge_shared_mutations,
    node_names,
  };

  let input = PruneInput {
    graph,
    alphabet,
    sequences,
    branch_lengths: branch_lengths_input,
  };

  progress.check_cancelled()?;
  progress.report("Pruning", 0.4, "");
  let output = pipeline::run(&params, input, &names)?;
  let pipeline::PruneOutput {
    graph,
    gtr,
    partitions,
    names,
    branch_lengths: branch_lengths_opt,
  } = output;
  // Share the sparse sequence partition Arc (no sequence data copied) and the fitted GTR into the
  // value-shaped result. The tree writers still read the partition and model from the graph data
  // slot until that read moves onto the result value; both copies leave the graph then.
  let seq = partitions.first().map(Arc::clone);
  let mut graph = graph.map_data(PruneGraphData::new(gtr.clone(), partitions));
  let topology_order = args
    .topology_order
    .resolve_topology_order(&graph, &names, Some(input_order))?;
  topology_order.apply(&mut graph, &names, &branch_lengths_opt)?;
  progress.report("Writing output", 0.8, "");

  // Gather the per-node name/confidence and per-edge branch length off the ordered tree into keyed
  // value maps the output writers consume. The name and branch length come from the post-topology
  // maps the pipeline returns (topology ordering only permutes keys, so their values still match the
  // pruned tree); the input branch support comes from the parse-time confidence map keyed by node,
  // which holds `None` for any node the pipeline created after the parse.
  let nodes: BTreeMap<GraphNodeKey, PruneNodeOut> = graph
    .get_nodes()
    .iter()
    .map(|node| {
      let node = node.read_arc();
      let key = node.key();
      let confidence = confidences.get(&key).copied().flatten();
      (
        key,
        PruneNodeOut {
          name: names[&key].clone(),
          confidence,
        },
      )
    })
    .collect();
  let edges: BTreeMap<GraphEdgeKey, EdgeOut> = graph
    .get_edges()
    .iter()
    .map(|edge| {
      let key = edge.read_arc().key();
      (
        key,
        EdgeOut {
          branch_length: branch_lengths_opt[&key],
        },
      )
    })
    .collect();
  let branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>> =
    edges.iter().map(|(key, edge)| (*key, edge.branch_length)).collect();

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::Gtr) {
    match graph.data().gtr.as_ref() {
      Some(gtr) => {
        let gtr_output = GtrOutput::new(gtr, GtrModelName::JC69);
        write_gtr_json(&gtr_output, path)?;
      },
      None if args.output_gtr.is_some() => {
        return make_error!(
          "GTR output requested but no GTR model was fitted. Provide sequence alignment input with --aln."
        );
      },
      None => warn!("Skipping GTR output: no GTR model was fitted (provide sequence alignment input with --aln)"),
    }
  }

  // Gather the per-node/per-edge sequence and mutation values off the partition into plain value maps
  // the output writers consume, taking the partition read out of the serialization path. Only the
  // auspice, phyloxml, and MAT writers read these maps; prune's default Newick/Nexus outputs carry no
  // mutation comments (empty comment provider), so they never read a mutation. Prune also applies its
  // final `--prune-empty` and `--merge-shared-mutations` topology edits without a following marginal
  // pass, leaving output-tree edges whose `subs_ml` was never populated. Gathering unconditionally
  // would read those unpopulated edges for outputs that never serialize them, so gather only when a
  // map-consuming tree output is requested.
  let maps = if resolved.tree_outputs.keys().any(prune_output_consumes_maps) {
    gather_prune_output_maps(&graph)?
  } else {
    PruneOutputMaps::default()
  };

  if !resolved.tree_outputs.is_empty() {
    write_prune_tree_outputs(
      &graph,
      &nodes,
      &branch_lengths,
      &maps,
      &resolved.tree_outputs,
      &CommentProviders::new(),
    )?;
  }

  progress.report("Done", 1.0, "");
  Ok(PruneResult {
    graph,
    nodes,
    edges,
    seq,
    gtr,
  })
}

/// Whether a tree-output kind reads the gathered prune value maps.
///
/// The auspice, phyloxml, and MAT writers build their nodes and edges from the value maps; the Newick,
/// Nexus, Graphviz, and internal-graph writers do not (prune supplies no mutation comment provider).
fn prune_output_consumes_maps(kind: &TreeWriteKind) -> bool {
  matches!(
    kind,
    TreeWriteKind::Auspice
      | TreeWriteKind::Phyloxml
      | TreeWriteKind::PhyloxmlJson
      | TreeWriteKind::MatPb
      | TreeWriteKind::MatJson
  )
}

/// Gather the per-node nucleotide sequences, root sequence, and per-edge nucleotide mutations the tree
/// writers read off the prune partition.
///
/// The per-edge mutation map is keyed by the inbound edge of every node reached on a walk from the
/// single root, i.e. exactly the output-tree edges the tree, MAT, and Newick-comment writers traverse.
/// Prune applies its final `--prune-empty` and `--merge-shared-mutations` topology edits without a
/// following marginal pass, so the node and edge stores can still hold detached nodes and orphan edges
/// whose `subs_ml` was never populated. Those never appear on the tree walk, so reading
/// `edge_mutations` only for reached edges avoids touching an unpopulated edge.
pub(crate) fn gather_prune_output_maps(graph: &GraphAncestral<PruneGraphData>) -> Result<PruneOutputMaps, Report> {
  let Some(partition) = graph.data().partitions.first() else {
    return Ok(PruneOutputMaps::default());
  };
  let partition = partition.read_arc();
  let root_sequence = Some(partition.root_sequence(graph)?);
  let node_sequences = graph
    .get_nodes()
    .iter()
    .map(|node| {
      let key = node.read_arc().key();
      (key, partition.node_sequence(key))
    })
    .collect();
  let mut edge_mutations = BTreeMap::new();
  let root = graph
    .get_exactly_one_root()
    .wrap_err("When gathering prune tree mutations")?;
  let mut queue = VecDeque::from([Arc::clone(&root)]);
  while let Some(node) = queue.pop_front() {
    for (child, edge) in graph.children_of(&node.read_arc()) {
      let edge_key = edge.read_arc().key();
      edge_mutations.insert(edge_key, partition.edge_mutations(graph, edge_key, MutationTrack::Nucleotide)?);
      queue.push_back(child);
    }
  }
  Ok(PruneOutputMaps {
    root_sequence,
    node_sequences,
    edge_mutations,
  })
}

fn validate_args(args: &TreetimePruneArgs) -> Result<(), Report> {
  if args.prune_empty && args.alignment.alignment.is_empty() {
    return make_error!(
      "The --prune-empty requires --aln. Without sequence data, it's not possible to determine which branches lack mutations."
    );
  }

  if args.merge_shared_mutations && args.alignment.alignment.is_empty() {
    return make_error!(
      "The --merge-shared-mutations requires --aln. Without sequence data, it's not possible to determine which branches share mutations."
    );
  }

  Ok(())
}

fn leaf_order<N, E, D>(
  graph: &Graph<N, E, D>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Result<Vec<String>, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Sync + Send,
{
  graph
    .get_leaves()
    .into_iter()
    .map(|leaf| {
      let key = leaf.read_arc().key();
      names[&key]
        .clone()
        .ok_or_else(|| crate::make_report!("Leaf node {key} has no name"))
    })
    .collect()
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
