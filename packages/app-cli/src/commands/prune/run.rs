use crate::commands::prune::args::TreetimePruneArgs;
use crate::commands::shared::resolve_outputs::ResolveOutputs;
use app_output::output_plan::OutputSelection;
use app_output::prune_result::{EdgeOut, PruneNodeOut, PruneOutputMaps, PruneResult};
use app_output::prune_tree_output::write_prune_tree_outputs;
use eyre::{Report, WrapErr};
use itertools::Itertools;
use log::warn;
use maplit::btreeset;
use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::path::PathBuf;
use treetime::alphabet::alphabet::Alphabet;
use treetime::ancestral::pipeline::SparseReconstruction;
use treetime::gtr::get_gtr::{GtrModelName, GtrOutput, write_gtr_json};
use treetime::make_error;
use treetime::prune::pipeline::{self, PruneInput, PruneParams};
use treetime::seq::mutation::MutationTrack;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::fasta::read_many_fasta_path;
use treetime_io::graph::TreeWriteKind;
use treetime_io::nwk::CommentProviders;
use treetime_io::nwk::nwk_read_file;
use treetime_io::parse_delimited::{parse_delimited_file, parse_delimited_str};
use treetime_primitives::AlignmentRecord;

pub fn run_prune(
  args: &TreetimePruneArgs,
  cancel: &dyn treetime::cancel::Cancel,
  progress: &dyn treetime::progress::ProgressSink,
) -> Result<PruneResult, Report> {
  validate_args(args)?;

  cancel.check()?;
  progress.report("Reading input", 0.0, "");

  let parse = nwk_read_file(args.tree())?;
  let confidences = parse.confidences();
  let names = parse.names();
  let graph: Graph = parse.graph;
  let branch_lengths_input = parse.branch_lengths;
  let input_order = leaf_order(&graph, &names)?;
  let alphabet = Alphabet::new(args.alphabet_args.alphabet_name().unwrap_or_default())?;

  let resolved = args.resolve_outputs()?;

  let needs_sequences = args.prune_empty || args.merge_shared_mutations;
  let sequences = if needs_sequences && !args.alignment.alignment.is_empty() {
    let records = read_many_fasta_path(&args.alignment.alignment, &alphabet)?;
    Some(records.into_iter().map(AlignmentRecord::from).collect())
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

  // Gather the per-node/per-edge sequence and mutation values off the pipeline-local partition into
  // plain value maps the output writers consume, taking the partition read out of the serialization
  // path. Only the auspice and MAT writers read these maps; prune's default Newick/Nexus
  // outputs carry no mutation comments (empty comment provider), so they never read a mutation. Prune
  // also applies its final `--prune-empty` and `--merge-shared-mutations` topology edits without a
  // following marginal pass, leaving output-tree edges whose `subs_ml` was never populated. Gathering
  // unconditionally would read those unpopulated edges for outputs that never serialize them, so
  // gather only when a map-consuming tree output is requested. Node and edge keys stay stable through
  // topology ordering, so gathering before it is bit-identical.
  let maps = if resolved.tree_outputs.keys().any(prune_output_consumes_maps) {
    gather_prune_output_maps(&graph, &partitions)?
  } else {
    PruneOutputMaps::default()
  };

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
      &maps,
      &resolved.tree_outputs,
      &CommentProviders::new(),
    )?;
  }

  progress.report("Done", 1.0, "");
  Ok(PruneResult { graph, nodes, edges })
}

/// Whether a tree-output kind reads the gathered prune value maps.
///
/// The auspice and MAT writers build their nodes and edges from the value maps; the Newick,
/// Nexus, Graphviz, and internal-graph writers do not (prune supplies no mutation comment provider).
fn prune_output_consumes_maps(kind: &TreeWriteKind) -> bool {
  matches!(
    kind,
    TreeWriteKind::Auspice | TreeWriteKind::MatPb | TreeWriteKind::MatJson
  )
}

#[allow(
  clippy::expect_used,
  reason = "expect on a value an upstream invariant guarantees is present"
)]
/// Gather the root sequence and per-edge nucleotide mutations the tree writers read off the prune
/// partition.
///
/// The per-edge mutation map is keyed by the inbound edge of every node reached on a walk from the
/// single root, i.e. exactly the output-tree edges the tree, MAT, and Newick-comment writers traverse.
/// Prune applies its final `--prune-empty` and `--merge-shared-mutations` topology edits without a
/// following marginal pass, so the node and edge stores can still hold detached nodes and orphan edges
/// whose `subs_ml` was never populated. Those never appear on the tree walk, so reading
/// `edge_mutations` only for reached edges avoids touching an unpopulated edge.
pub(crate) fn gather_prune_output_maps(
  graph: &Graph,
  partitions: &[SparseReconstruction],
) -> Result<PruneOutputMaps, Report> {
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

fn leaf_order(graph: &Graph, names: &BTreeMap<GraphNodeKey, Option<String>>) -> Result<Vec<String>, Report> {
  graph
    .get_leaves()
    .map(|leaf| {
      let key = leaf.key();
      names[&key]
        .clone()
        .ok_or_else(|| treetime::make_report!("Leaf node {key} has no name"))
    })
    .collect()
}

#[allow(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
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
