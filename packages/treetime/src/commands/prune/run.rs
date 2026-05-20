use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::fitch::create_fitch_partition;
use crate::commands::prune::args::TreetimePruneArgs;
use crate::gtr::get_gtr::{GtrModelName, get_gtr_by_name, log_gtr, write_gtr_json};
use crate::make_error;
use crate::optimize::topology::merge_shared_mutations::merge_shared_mutation_branches;
use crate::partition::marginal_sparse::PartitionMarginalSparse;
use crate::payload::ancestral::GraphAncestral;
use crate::prune::prune::prune_nodes;
use eyre::Report;
use itertools::Itertools;
use maplit::btreeset;
use parking_lot::RwLock;
use std::collections::BTreeSet;
use std::path::PathBuf;
use std::sync::Arc;
use treetime_io::fasta::read_many_fasta;
use treetime_io::graph::write_graph_files;
use treetime_io::nwk::nwk_read_file;
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

    let fitch = create_fitch_partition(&graph, 0, alphabet, &aln)?;
    let gtr = get_gtr_by_name(GtrModelName::JC69)?;
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

  prune_nodes(&mut graph, &partitions, *prune_short, *prune_empty, &node_names)?;

  if *merge_shared_mutations {
    merge_shared_mutation_branches(&mut graph, &partitions)?;
    graph.build()?;
  }

  write_graph_files(outdir, "pruned_tree", &graph)?;

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
