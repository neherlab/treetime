use crate::commands::prune::args::TreetimePruneArgs;
use crate::commands::shared::alignment::read_alignment;
use crate::commands::shared::leaf_order::leaf_order;
use crate::commands::shared::resolve_outputs::ResolveOutputs;
use app_output::annotated_graph::{AnnotatedGraph, Divergence, TreeSequences};
use app_output::mutation_filter::UnknownMutationFilter;
use app_output::output_plan::{CommandKind, OutputSelection, ResolvedOutputs, TreeWriteKind};
use app_output::tree_output::{tree_view_for_outputs, write_graph_outputs, write_tree_outputs};
use eyre::Report;
use maplit::btreeset;
use std::collections::{BTreeMap, BTreeSet};
use std::path::PathBuf;
use treetime::alphabet::alphabet::Alphabet;
use treetime::cancel::Cancel;
use treetime::gtr::get_gtr::{GtrModelName, GtrOutput};
use treetime::partition::marginal::sparse::partition::PartitionMarginalSparse;
use treetime::progress::{LogSink, StageSink};
use treetime::progress_warn;
use treetime::prune::pipeline::{self, PruneInput, PruneParams};
use treetime::seq::mutation::{Mutation, MutationTrack};
use treetime::{make_error, make_report};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::name_list::{name_list_read_file, name_list_read_str};
use treetime_io::nwk::{NwkStyle, nwk_read_file};
use treetime_primitives::{AlignmentRecord, Seq};
use treetime_utils::io::json::{JsonPretty, json_write_file};

pub fn run_prune(
  args: &TreetimePruneArgs,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<(), Report> {
  validate_args(args)?;

  cancel.check()?;
  stages.report("Reading input", 0.0, "");

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
    let records = read_alignment(&args.alignment.alignment, &alphabet)?;
    Some(records.into_iter().map(AlignmentRecord::from).collect())
  } else {
    None
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
    node_names,
  };

  let unknown = alphabet.unknown();
  let input = PruneInput {
    graph,
    alphabet,
    sequences,
    branch_lengths: branch_lengths_input,
  };

  cancel.check()?;
  stages.report("Pruning", 0.4, "");
  let output = pipeline::run(&params, input, &names, cancel, log).map_err(|err| err.into_report())?;
  let pipeline::PruneOutput {
    mut graph,
    gtr,
    partitions,
    names,
    branch_lengths: branch_lengths_opt,
  } = output;

  let maps = if resolved.tree_outputs.keys().copied().any(prune_output_consumes_maps) {
    gather_prune_output_maps(&graph, &partitions, UnknownMutationFilter::hiding_unknown(unknown))?
  } else {
    PruneOutputMaps::default()
  };

  let topology_order = args
    .topology_order
    .resolve_topology_order(&graph, &names, Some(input_order))?;
  topology_order.apply(&mut graph, &names, &branch_lengths_opt)?;
  stages.report("Writing output", 0.8, "");

  let branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>> = graph
    .get_edges()
    .map(|edge| (edge.key(), branch_lengths_opt[&edge.key()]))
    .collect();

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::Gtr) {
    match gtr.as_ref() {
      Some(gtr) => {
        let gtr_output = GtrOutput::builder().gtr(gtr).model_name(GtrModelName::JC69).build();
        json_write_file(path, &gtr_output, JsonPretty(true))?;
      },
      None if args.output_gtr.is_some() => {
        return make_error!(
          "GTR output requested but no GTR model was fitted. Prune fits the model from the alignment only for --prune-empty or --merge-shared-mutations."
        );
      },
      None => progress_warn!(
        log,
        "Skipping GTR output: no GTR model was fitted (prune fits it from the alignment only for --prune-empty or --merge-shared-mutations)"
      ),
    }
  }

  write_prune_trees(&graph, &names, &branch_lengths, &confidences, &maps, &resolved, log)?;

  stages.report("Done", 1.0, "");
  Ok(())
}

fn write_prune_trees(
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  branch_support: &BTreeMap<GraphNodeKey, Option<f64>>,
  maps: &PruneOutputMaps,
  resolved: &ResolvedOutputs,
  log: &dyn LogSink,
) -> Result<(), Report> {
  let annotated = AnnotatedGraph {
    graph,
    names,
    divergence_branch_lengths: branch_lengths,
    time_branch_lengths: None,
    divergence: Divergence::CumulativeBranchLength,
    branch_support: Some(branch_support),
    sequences: maps.root_sequence.as_ref().map(|root_sequence| TreeSequences {
      root_sequence,
      edge_mutations: &maps.edge_mutations,
      mutation_counts: None,
      amino_acids: None,
    }),
    dates: None,
    traits: None,
  };
  write_graph_outputs(&annotated, &resolved.tree_outputs)?;
  if let Some(tree) = tree_view_for_outputs(&annotated, resolved)? {
    write_tree_outputs(&tree, &resolved.tree_outputs, CommandKind::Prune, log)?;
  }
  Ok(())
}

fn prune_output_consumes_maps(kind: TreeWriteKind) -> bool {
  match kind {
    TreeWriteKind::Auspice | TreeWriteKind::MatPb | TreeWriteKind::MatJson => true,
    TreeWriteKind::Nwk(style) | TreeWriteKind::Nexus(style) => style != NwkStyle::Plain,
    TreeWriteKind::GraphJson | TreeWriteKind::Dot => false,
  }
}

#[derive(Debug, Default)]
struct PruneOutputMaps {
  root_sequence: Option<Seq>,
  edge_mutations: BTreeMap<GraphEdgeKey, Vec<Mutation>>,
}

fn gather_prune_output_maps(
  graph: &Graph,
  partitions: &[PartitionMarginalSparse],
  filter: UnknownMutationFilter,
) -> Result<PruneOutputMaps, Report> {
  let Some(partition) = partitions.first() else {
    return Ok(PruneOutputMaps::default());
  };
  let edge_mutations = graph
    .get_edges()
    .map(|edge| {
      let edge_key = edge.key();
      Ok((
        edge_key,
        partition.edge_fitch_mutations(edge_key, &MutationTrack::Nucleotide)?,
      ))
    })
    .collect::<Result<BTreeMap<_, _>, Report>>()?;
  Ok(PruneOutputMaps {
    root_sequence: Some(partition.fitch_root_sequence()),
    edge_mutations: filter.reported_edge_mutations(graph, edge_mutations)?,
  })
}

fn validate_args(args: &TreetimePruneArgs) -> Result<(), Report> {
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

fn ascii_delimiter(delimiter: char) -> Result<u8, Report> {
  u8::try_from(delimiter)
    .ok()
    .filter(u8::is_ascii)
    .ok_or_else(|| make_report!("the delimiter must be an ASCII character, got '{delimiter}'"))
}
