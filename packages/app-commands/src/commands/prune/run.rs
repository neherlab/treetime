use crate::commands::prune::args::TreetimePruneArgs;
use crate::commands::shared::leaf_order::leaf_order;
use crate::commands::shared::partition::{log_partition_layouts, resolve_alignment_layouts};
use crate::commands::shared::tree_input::read_input_tree;
use app_output::annotated_graph::{AnnotatedGraph, Divergence, TreeSequences};
use app_output::command_outputs::CommandKind;
use app_output::mutation_filter::UnknownMutationFilter;
use app_output::output_plan::{ResolvedOutputs, TreeWriteKind};
use app_output::tree_output::{write_run_tree_outputs, write_tree_outputs_per_alignment};
use eyre::Report;
use itertools::Itertools;
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;
use treetime::cancel::Cancel;
use treetime::partition::result::AlignmentSequences;
use treetime::progress::{LogSink, StageSink};
use treetime::prune::pipeline::{self, PruneInput, PruneParams};
use treetime::seq::mutation::Mutation;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::EdgeMap;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::node_names::NodeNames;
use treetime_graph::pair_by_name::key_by_name;
use treetime_io::name_list::name_list_read_file;
use treetime_primitives::Seq;

pub fn run_prune(
  args: &TreetimePruneArgs,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<(), Report> {
  let resolved = args.outputs.resolve()?;

  cancel.check()?;
  stages.report("Reading input", 0.0, "");

  let parse = read_input_tree(&args.tree, args.newick_dialect, log)?;
  let names = parse.names;
  let graph: Graph = parse.graph;
  let branch_lengths_input = parse.branch_lengths;
  let input_order = leaf_order(&graph);
  let sources = args
    .sequences
    .read(&graph, &names, log)?
    .into_iter()
    .map(|source| source.source)
    .collect_vec();
  let layouts = resolve_alignment_layouts(&sources.iter().map(|source| source.length).collect_vec())?;
  log_partition_layouts(&layouts, log);

  let node_names = node_names_to_remove(&args.remove, args.remove_file.as_deref())?;

  let params = PruneParams {
    prune_short: args.collapse_short,
    prune_empty: args.collapse_empty,
    merge_shared_mutations: args.merge_shared_mutations,
    node_keys: resolve_node_names(&node_names, &graph, &names),
  };

  let input = PruneInput {
    graph,
    names,
    alignments: sources,
    layouts,
    branch_lengths: branch_lengths_input,
  };

  cancel.check()?;
  stages.report("Pruning", 0.4, "");
  let output = pipeline::run(&params, input, cancel, log).map_err(|err| err.into_report())?;
  let pipeline::PruneOutput {
    mut graph,
    alignments: mut results,
    names,
    branch_lengths: branch_lengths_opt,
    ..
  } = output;
  let maps = if resolved.tree_write_kinds().any(TreeWriteKind::is_per_alignment) {
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

  let topology_order = args
    .topology_order
    .resolve(&graph, &names, Some(input_order), args.newick_dialect)?;
  topology_order.apply(&mut graph, &names, &branch_lengths_opt)?;
  stages.report("Writing output", 0.8, "");

  let branch_lengths: EdgeMap<Option<f64>> = graph
    .get_edges()
    .map(|edge| (edge.key(), branch_lengths_opt[&edge.key()]))
    .collect();

  let sequences = results
    .iter()
    .filter_map(|result| {
      let maps = maps.get(&result.number)?;
      let sequences = TreeSequences {
        alphabet: &result.alphabet,
        root_sequence: &maps.root_sequence,
        edge_mutations: &maps.edge_mutations,
        amino_acids: None,
        comment_positions: args.positions,
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
  names: &NodeNames,
  branch_lengths: &EdgeMap<Option<f64>>,
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
  write_tree_outputs_per_alignment(&annotated, sequences, resolved, CommandKind::Prune, log)
}

fn prune_output_maps(
  graph: &Graph,
  sequences: AlignmentSequences,
  filter: UnknownMutationFilter,
) -> Result<PruneOutputMaps, Report> {
  Ok(PruneOutputMaps {
    root_sequence: sequences.root_sequence,
    edge_mutations: filter.reported_edge_mutations(graph, sequences.edge_mutations)?,
  })
}

struct PruneOutputMaps {
  root_sequence: Seq,
  edge_mutations: EdgeMap<Vec<Mutation>>,
}

fn node_names_to_remove(remove: &[String], remove_file: Option<&Path>) -> Result<BTreeSet<String>, Report> {
  let from_file = remove_file.map(name_list_read_file).transpose()?;
  Ok(remove.iter().cloned().chain(from_file.into_iter().flatten()).collect())
}

pub(crate) fn resolve_node_names(
  node_names: &BTreeSet<String>,
  graph: &Graph,
  names: &NodeNames,
) -> BTreeSet<GraphNodeKey> {
  let node_of = key_by_name(graph.get_nodes().map(|node| node.key()), names);
  node_names
    .iter()
    .filter_map(|name| node_of.get(name.as_str()).copied())
    .collect()
}
