use crate::alphabet::alphabet::Alphabet;
use crate::commands::prune::args::TreetimePruneArgs;
use crate::commands::prune::result::{EdgeOut, PruneGraphData, PruneNodeOut, PruneResult};
use crate::commands::shared::output::OutputSelection;
use crate::commands::shared::resolve_outputs::ResolveOutputs;
use crate::commands::shared::tree_output::write_prune_tree_outputs;
use crate::gtr::get_gtr::{GtrModelName, GtrOutput, write_gtr_json};
use crate::make_error;
use crate::prune::pipeline::{self, PruneInput, PruneParams};
use eyre::Report;
use itertools::Itertools;
use log::warn;
use maplit::btreeset;
use std::collections::{BTreeMap, BTreeSet};
use std::path::PathBuf;
use std::sync::Arc;
use treetime_graph::edge::{GraphEdge, GraphEdgeKey};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, GraphNodeKey, Named};
use treetime_graph::value_maps::node_names;
use treetime_io::fasta::read_many_fasta;
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

  let graph: GraphAncestral = nwk_read_file(args.tree())?;
  let input_order = leaf_order(&graph)?;
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
  };

  progress.check_cancelled()?;
  progress.report("Pruning", 0.4, "");
  let output = pipeline::run(&params, input)?;
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
  let topology_order = args.topology_order.resolve_topology_order(&graph, Some(input_order))?;
  topology_order.apply(&mut graph)?;
  progress.report("Writing output", 0.8, "");

  // Gather the per-node name/confidence and per-edge branch length off the ordered tree into keyed
  // value maps the output writers consume. The name and branch length come from the post-topology
  // maps the pipeline returns (topology ordering only permutes keys, so their values still match the
  // pruned tree); the input branch support is still read from the payload until it moves onto a map.
  let nodes: BTreeMap<GraphNodeKey, PruneNodeOut> = graph
    .get_nodes()
    .iter()
    .map(|node| {
      let node = node.read_arc();
      let key = node.key();
      let confidence = node.payload().read_arc().confidence;
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

  if !resolved.tree_outputs.is_empty() {
    write_prune_tree_outputs(
      &graph,
      &nodes,
      &branch_lengths,
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

fn leaf_order<N, E, D>(graph: &Graph<N, E, D>) -> Result<Vec<String>, Report>
where
  N: GraphNode + Named,
  E: GraphEdge,
  D: Sync + Send,
{
  let names = node_names(graph);
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
