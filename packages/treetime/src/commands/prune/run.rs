use crate::alphabet::alphabet::Alphabet;
use crate::commands::prune::args::TreetimePruneArgs;
use crate::commands::prune::result::PruneResult;
use crate::commands::shared::output::{CommandKind, OutputSelection};
use crate::gtr::get_gtr::{GtrModelName, GtrOutput, write_gtr_json};
use crate::make_error;
use crate::prune::pipeline::{self, PruneInput, PruneParams};
use eyre::Report;
use itertools::Itertools;
use maplit::btreeset;
use std::collections::BTreeSet;
use std::path::PathBuf;
use treetime_graph::edge::GraphEdge;
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, Named};
use treetime_io::fasta::read_many_fasta;
use treetime_io::graph::write_tree_outputs;
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

  let graph: GraphAncestral = nwk_read_file(&args.tree)?;
  let input_order = leaf_order(&graph)?;
  let alphabet = Alphabet::new(args.alphabet_args.alphabet.unwrap_or_default())?;

  let resolved = args.output.resolve(
    CommandKind::Prune,
    &graph,
    &[(OutputSelection::Gtr, args.output_gtr.as_deref())],
    Some(input_order),
  )?;

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

  progress.report("Writing output", 0.8, "");

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::Gtr) {
    let gtr = output.gtr.as_ref().ok_or_else(|| {
      crate::make_report!(
        "GTR output requested but no GTR model was fitted. Provide sequence alignment input with --aln."
      )
    })?;
    let gtr_output = GtrOutput::new(gtr, GtrModelName::JC69);
    write_gtr_json(&gtr_output, path)?;
  }

  if !resolved.tree_outputs.is_empty() {
    let plan = resolved.topology_order.plan(&output.graph)?;
    let ordered = plan.ordered_graph(&output.graph)?;
    write_tree_outputs(&ordered, &resolved.tree_outputs, &CommentProviders::new(), None)?;
  }

  progress.report("Done", 1.0, "");
  Ok(PruneResult { graph: output.graph })
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
  graph
    .get_leaves()
    .into_iter()
    .map(|leaf| {
      let leaf = leaf.read_arc();
      leaf
        .payload()
        .read_arc()
        .name()
        .map(|name| name.as_ref().to_owned())
        .ok_or_else(|| crate::make_report!("Leaf node {} has no name", leaf.key()))
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
