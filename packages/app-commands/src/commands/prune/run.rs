use crate::commands::prune::args::TreetimePruneArgs;
use crate::commands::shared::gtr_output::{partition_models, write_gtr_output};
use crate::commands::shared::leaf_order::leaf_order;
use crate::commands::shared::output_args::ResolveOutputs;
use crate::commands::shared::sequence_inputs::{SourceSettings, read_sequence_sources};
use crate::commands::shared::tree_input::read_input_tree;
use app_output::annotated_graph::{AnnotatedGraph, Divergence, TreeSequences};
use app_output::mutation_filter::UnknownMutationFilter;
use app_output::output_plan::{CommandKind, ResolvedOutputs, TreeWriteKind};
use app_output::tree_output::{write_run_tree_outputs, write_tree_outputs_per_partition};
use eyre::Report;
use maplit::btreeset;
use std::collections::{BTreeMap, BTreeSet};
use std::path::PathBuf;
use treetime::alphabet::alphabet::Alphabet;
use treetime::ancestral::attach::MissingLeafLimit;
use treetime::cancel::Cancel;
use treetime::partition::result::PartitionSequences;
use treetime::progress::{LogSink, StageSink};
use treetime::prune::pipeline::{self, PruneInput, PruneParams};
use treetime::seq::mutation::Mutation;
use treetime::{make_error, make_report};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::pair_by_name::first_key_by_name;
use treetime_io::name_list::{name_list_read_file, name_list_read_str};
use treetime_primitives::Seq;

pub fn run_prune(
  args: &TreetimePruneArgs,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<(), Report> {
  validate_args(args)?;

  cancel.check()?;
  stages.report("Reading input", 0.0, "");

  let parse = read_input_tree(args.tree(), args.tree_dialect.dialect(), log)?;
  let names = parse.names();
  let graph: Graph = parse.graph;
  let branch_lengths_input = parse.branch_lengths;
  let input_order = leaf_order(&graph);
  let alphabet = Alphabet::new(args.alphabet_args.alphabet_name().unwrap_or_default())?;

  let resolved = args.resolve_outputs()?;

  let needs_sequences = args.prune_empty || args.merge_shared_mutations;
  let sources = if needs_sequences && !args.alignment.alignment.is_empty() {
    let settings = SourceSettings {
      alphabet: &alphabet,
      gap_fill: None,
      missing_leaves: MissingLeafLimit::Enforced,
    };
    read_sequence_sources(&args.alignment, &settings, &graph, &names, log)?
      .into_iter()
      .map(|source| source.source)
      .collect()
  } else {
    vec![]
  };

  let node_names = parse_node_names(
    args.prune_nodes_list.as_deref(),
    args.prune_nodes_list_delimiter,
    args.prune_nodes_list_file.as_ref(),
    args.prune_nodes_list_file_delimiter,
  )?;

  let params = PruneParams {
    prune_short: args.prune_short,
    prune_empty: args.prune_empty,
    merge_shared_mutations: args.merge_shared_mutations,
    node_keys: resolve_node_names(&node_names, &graph, &names),
  };

  let input = PruneInput {
    graph,
    sources,
    branch_lengths: branch_lengths_input,
  };

  cancel.check()?;
  stages.report("Pruning", 0.4, "");
  let output = pipeline::run(&params, input, &names, cancel, log).map_err(|err| err.into_report())?;
  let pipeline::PruneOutput {
    mut graph,
    mut results,
    names,
    branch_lengths: branch_lengths_opt,
  } = output;
  let models = partition_models(&results, None);
  let maps = if resolved.tree_write_kinds().any(TreeWriteKind::is_per_partition) {
    results
      .iter_mut()
      .filter_map(|result| {
        let sequences = result.sequences.take()?;
        let filter = UnknownMutationFilter::hiding_unknown(result.alphabet.unknown());
        Some(prune_output_maps(&graph, sequences, filter).map(|maps| (result.number, maps)))
      })
      .collect::<Result<BTreeMap<_, _>, Report>>()?
  } else {
    BTreeMap::new()
  };

  let topology_order =
    args
      .topology_order
      .resolve_topology_order(&graph, &names, Some(input_order), args.tree_dialect.dialect())?;
  topology_order.apply(&mut graph, &names, &branch_lengths_opt)?;
  stages.report("Writing output", 0.8, "");

  let branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>> = graph
    .get_edges()
    .map(|edge| (edge.key(), branch_lengths_opt[&edge.key()]))
    .collect();

  write_gtr_output(
    &resolved,
    models,
    "no GTR model was fitted (prune fits it from the alignment only for --prune-empty or --merge-shared-mutations)",
    log,
  )?;

  let sequences = results
    .iter()
    .filter_map(|result| {
      let maps = maps.get(&result.number)?;
      let sequences = TreeSequences {
        alphabet: &result.alphabet,
        root_sequence: &maps.root_sequence,
        edge_mutations: &maps.edge_mutations,
        amino_acids: None,
      };
      Some((result.number, sequences))
    })
    .collect();
  write_prune_trees(&graph, &names, &branch_lengths, &sequences, &resolved, log)?;

  stages.report("Done", 1.0, "");
  Ok(())
}

fn write_prune_trees(
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  sequences: &BTreeMap<usize, TreeSequences<'_>>,
  resolved: &ResolvedOutputs,
  log: &dyn LogSink,
) -> Result<(), Report> {
  let annotated = AnnotatedGraph {
    graph,
    names,
    divergence_branch_lengths: branch_lengths,
    time_branch_lengths: None,
    divergence: Divergence::CumulativeBranchLength,
    sequences: None,
    dates: None,
    traits: None,
  };
  write_run_tree_outputs(&annotated, resolved, CommandKind::Prune, log)?;
  write_tree_outputs_per_partition(&annotated, sequences, resolved, CommandKind::Prune, log)
}

fn prune_output_maps(
  graph: &Graph,
  sequences: PartitionSequences,
  filter: UnknownMutationFilter,
) -> Result<PruneOutputMaps, Report> {
  Ok(PruneOutputMaps {
    root_sequence: sequences.root_sequence,
    edge_mutations: filter.reported_edge_mutations(graph, sequences.edge_mutations)?,
  })
}

struct PruneOutputMaps {
  root_sequence: Seq,
  edge_mutations: BTreeMap<GraphEdgeKey, Vec<Mutation>>,
}

fn validate_args(args: &TreetimePruneArgs) -> Result<(), Report> {
  args.alignment.paths()?;

  if args.prune_empty && args.alignment.alignment.is_empty() {
    return make_error!(
      "The --prune-empty requires --alignment. Without sequence data, it's not possible to determine which branches lack mutations."
    );
  }

  if args.merge_shared_mutations && args.alignment.alignment.is_empty() {
    return make_error!(
      "The --merge-shared-mutations requires --alignment. Without sequence data, it's not possible to determine which branches share mutations."
    );
  }

  Ok(())
}

fn parse_node_names(
  prune_nodes_list: Option<&str>,
  prune_nodes_list_delimiter: char,
  prune_nodes_list_file: Option<&PathBuf>,
  prune_nodes_list_file_delimiter: char,
) -> Result<BTreeSet<String>, Report> {
  let mut node_names = btreeset! {};

  if let Some(prune_nodes_list) = prune_nodes_list {
    node_names.extend(name_list_read_str(
      prune_nodes_list,
      ascii_delimiter(prune_nodes_list_delimiter)?,
    )?);
  }

  if let Some(prune_nodes_list_file) = prune_nodes_list_file {
    node_names.extend(name_list_read_file(
      prune_nodes_list_file,
      ascii_delimiter(prune_nodes_list_file_delimiter)?,
    )?);
  }

  Ok(node_names)
}

pub(crate) fn resolve_node_names(
  node_names: &BTreeSet<String>,
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> BTreeSet<GraphNodeKey> {
  let first_node_of = first_key_by_name(graph.get_nodes().map(|node| node.key()), names);
  node_names
    .iter()
    .filter_map(|name| first_node_of.get(name.as_str()).copied())
    .collect()
}

fn ascii_delimiter(delimiter: char) -> Result<u8, Report> {
  u8::try_from(delimiter)
    .ok()
    .filter(u8::is_ascii)
    .ok_or_else(|| make_report!("the delimiter must be an ASCII character, got '{delimiter}'"))
}
