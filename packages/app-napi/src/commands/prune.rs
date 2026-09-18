//! N-API `prune` request shape and orchestration.

use crate::commands::support::{default_output_plan, default_topology_order};
use app_output::output_plan::CommandKind;
use app_output::output_plan::OutputSelection;
use app_output::prune_result::{EdgeOut, PruneNodeOut, PruneOutputMaps, PruneResult};
use app_output::prune_tree_output::write_prune_tree_outputs;
use eyre::{Report, WrapErr};
use itertools::Itertools;
use maplit::btreeset;
use serde::Deserialize;
use smart_default::SmartDefault;
use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::path::Path;
use treetime::alphabet::alphabet::{Alphabet, AlphabetName};
use treetime::ancestral::pipeline::SparseReconstruction;
use treetime::cancel::Cancel;
use treetime::gtr::get_gtr::{GtrModelName, GtrOutput, write_gtr_json};
use treetime::make_error;
use treetime::progress::ProgressSink;
use treetime::prune::pipeline::{self, PruneInput, PruneParams};
use treetime::seq::mutation::MutationTrack;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::fasta::read_many_fasta_path;
use treetime_io::graph::TreeWriteKind;
use treetime_io::nwk::{CommentProviders, nwk_read_file};
use treetime_io::parse_delimited::{parse_delimited_file, parse_delimited_str};
use treetime_primitives::AlignmentRecord;

/// Tree-pruning request (openapi subset).
#[derive(Debug, SmartDefault, Deserialize)]
#[serde(default)]
pub struct PruneArgs {
  pub input_fastas: Vec<String>,
  pub tree: String,
  pub alphabet: Option<AlphabetName>,
  pub outdir: String,
  pub prune_short: Option<f64>,
  pub prune_empty: bool,
  pub merge_shared_mutations: bool,
  pub prune_nodes_list: Option<String>,
  #[default = ',']
  pub prune_nodes_list_delimiter: char,
  pub prune_nodes_list_file: Option<String>,
  #[default = '\n']
  pub prune_nodes_list_file_delimiter: char,
}

pub fn run_prune(args: &PruneArgs, cancel: &dyn Cancel, progress: &dyn ProgressSink) -> Result<PruneResult, Report> {
  validate_args(args)?;

  cancel.check()?;
  progress.report("Reading input", 0.0, "");

  let parse = nwk_read_file(Path::new(&args.tree))?;
  let confidences = parse.confidences();
  let names = parse.names();
  let graph: Graph = parse.graph;
  let branch_lengths_input = parse.branch_lengths;
  let alphabet = Alphabet::new(args.alphabet.unwrap_or_default())?;

  let resolved = default_output_plan(CommandKind::Prune, Path::new(&args.outdir))?;

  let needs_sequences = args.prune_empty || args.merge_shared_mutations;
  let sequences = if needs_sequences && !args.input_fastas.is_empty() {
    let paths: Vec<std::path::PathBuf> = args.input_fastas.iter().map(std::path::PathBuf::from).collect();
    let records = read_many_fasta_path(&paths, &alphabet)?;
    Some(records.into_iter().map(AlignmentRecord::from).collect())
  } else {
    None
  };

  let node_names = parse_node_names(
    args.prune_nodes_list.as_ref(),
    args.prune_nodes_list_delimiter,
    args.prune_nodes_list_file.as_deref(),
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

  cancel.check()?;
  progress.report("Pruning", 0.4, "");
  let output = pipeline::run(&params, input, &names, cancel).map_err(|err| err.into_report())?;
  let pipeline::PruneOutput {
    mut graph,
    gtr,
    partitions,
    names,
    branch_lengths: branch_lengths_opt,
  } = output;

  let maps = if resolved.tree_outputs.keys().any(prune_output_consumes_maps) {
    gather_prune_output_maps(&graph, &partitions)?
  } else {
    PruneOutputMaps::default()
  };

  default_topology_order().apply(&mut graph, &names, &branch_lengths_opt)?;
  progress.report("Writing output", 0.8, "");

  let nodes: BTreeMap<GraphNodeKey, PruneNodeOut> = graph
    .get_nodes()
    .map(|node| {
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
    .map(|edge| {
      let key = edge.key();
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
    match gtr.as_ref() {
      Some(gtr) => {
        let gtr_output = GtrOutput::new(gtr, GtrModelName::JC69);
        write_gtr_json(&gtr_output, path)?;
      },
      None => log::warn!("Skipping GTR output: no GTR model was fitted (provide sequence alignment input with --aln)"),
    }
  }

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
  Ok(PruneResult { graph, nodes, edges })
}

fn prune_output_consumes_maps(kind: &TreeWriteKind) -> bool {
  matches!(
    kind,
    TreeWriteKind::Auspice | TreeWriteKind::MatPb | TreeWriteKind::MatJson
  )
}

fn gather_prune_output_maps(graph: &Graph, partitions: &[SparseReconstruction]) -> Result<PruneOutputMaps, Report> {
  let Some(partition) = partitions.first() else {
    return Ok(PruneOutputMaps::default());
  };
  let root_sequence = Some(partition.root_sequence(graph)?);
  let mut edge_mutations = BTreeMap::new();
  let root_key = graph
    .get_exactly_one_root()
    .wrap_err("When gathering prune tree mutations")?
    .key();
  let mut queue = VecDeque::from([root_key]);
  while let Some(node_key) = queue.pop_front() {
    let node = graph.get_node(node_key).expect("Node from graph traversal must exist");
    for (child_key, edge_key) in graph.children_keys_of(node) {
      edge_mutations.insert(
        edge_key,
        partition.edge_mutations(edge_key, &MutationTrack::Nucleotide)?,
      );
      queue.push_back(child_key);
    }
  }
  Ok(PruneOutputMaps {
    root_sequence,
    edge_mutations,
  })
}

fn validate_args(args: &PruneArgs) -> Result<(), Report> {
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
  prune_nodes_list_file: Option<&str>,
  prune_nodes_list_file_delimiter: char,
) -> Result<BTreeSet<String>, Report> {
  let mut node_names = btreeset! {};

  if let Some(prune_nodes_list) = prune_nodes_list {
    let names: Vec<String> = parse_delimited_str(prune_nodes_list, prune_nodes_list_delimiter as u8).try_collect()?;
    node_names.extend(names);
  }

  if let Some(prune_nodes_list_file) = prune_nodes_list_file {
    let names: Vec<String> =
      parse_delimited_file(Path::new(prune_nodes_list_file), prune_nodes_list_file_delimiter as u8)?.try_collect()?;
    node_names.extend(names);
  }

  Ok(node_names)
}
