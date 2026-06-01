use crate::alphabet::alphabet::Alphabet;
use crate::commands::prune::args::TreetimePruneArgs;
use crate::commands::prune::result::PruneResult;
use crate::gtr::get_gtr::{GtrModelName, write_gtr_json};
use crate::make_error;
use crate::prune::pipeline::{self, PruneInput, PruneParams};
use eyre::Report;
use itertools::Itertools;
use maplit::btreeset;
use std::collections::BTreeSet;
use std::path::PathBuf;
use treetime_io::fasta::read_many_fasta;
use treetime_io::graph::write_graph_files;
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
  let alphabet = Alphabet::new(args.alphabet_args.alphabet.unwrap_or_default())?;

  let sequences = if !args.alignment.alignment.is_empty() {
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
  let outdir = &args.output.outdir;
  if let Some(gtr) = &output.gtr {
    write_gtr_json(gtr, GtrModelName::JC69, outdir, None)?;
  }
  write_graph_files(outdir, "pruned_tree", &output.graph)?;

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
