use crate::commands::optimize::args::TreetimeOptimizeArgs;
use crate::commands::shared::gtr_output::{partition_models, write_gtr_output};
use crate::commands::shared::output_args::DivergenceUnits;
use crate::commands::shared::output_args::ResolveOutputs;
use crate::commands::shared::sequence_inputs::{SourceSettings, read_sequence_sources, single_source};
use crate::commands::shared::tree_input::read_input_tree;
use app_output::annotated_graph::{AnnotatedGraph, Divergence, TreeSequences};
use app_output::augur_node_data_refine::{RefineRun, write_augur_node_data_refine};
use app_output::mutation_filter::UnknownMutationFilter;
use app_output::output_plan::{CommandKind, OutputSelection, ResolvedOutputs};
use app_output::tree_output::{tree_view_for_outputs, write_graph_outputs, write_tree_outputs};
use eyre::Report;
use itertools::Itertools;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use treetime::alphabet::alphabet::Alphabet;
use treetime::ancestral::attach::MissingLeafLimit;
use treetime::cancel::Cancel;
use treetime::optimize::params::OptimizeLoopSettings;
use treetime::optimize::pipeline::{self, OptimizeInput, OptimizeParams};
use treetime::partition::result::PartitionSequences;
use treetime::progress::{LogSink, StageSink};
use treetime::progress_info;
use treetime::seq::mutation::Mutation;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::Seq;
use treetime_utils::make_internal_error;

pub fn run_optimize(
  args: &TreetimeOptimizeArgs,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<(), Report> {
  args.alignment.required_paths()?;
  cancel.check()?;
  stages.report("Reading input", 0.0, "");

  let alphabet = Alphabet::new(args.alphabet_args.alphabet_name().unwrap_or_default())?;
  let nwk_parsed = read_input_tree(args.tree(), args.tree_dialect.dialect(), log)?;
  let names = nwk_parsed.names();
  let graph = nwk_parsed.graph;
  let branch_lengths = nwk_parsed.branch_lengths;
  let settings = SourceSettings {
    alphabet: &alphabet,
    gap_fill: Some(args.gap_fill_args.effective_gap_fill()),
    missing_leaves: MissingLeafLimit::Enforced,
  };
  let sources = read_sequence_sources(&args.alignment, &settings, &graph, &names, log)?;
  let source = single_source(sources)?.source;

  let resolved = args.resolve_outputs()?;

  let params = OptimizeParams {
    model: args.model_args.model_name(),
    dense: args.dense,
    initial_guess: args.branch_length_initial_guess,
    reroot_spec: args.reroot_spec(&graph, &names)?,
    loop_settings: OptimizeLoopSettings {
      max_iter: args.max_iter,
      dp: args.dp,
      damping: args.damping,
      opt_method: args.opt_method,
      indel_policy: args.indel_policy,
      topology_ops: args.topology_ops,
    },
  };

  let unknown = alphabet.unknown();
  let input = OptimizeInput {
    graph,
    source,
    branch_lengths,
  };

  let output = pipeline::run(&params, input, &names, cancel, stages, log).map_err(|err| err.into_report())?;
  let pipeline::OptimizeOutput {
    mut graph,
    results,
    branch_lengths,
    names,
  } = output;
  let models = partition_models(&results, None);
  let Ok(mut result) = results.into_iter().exactly_one() else {
    return make_internal_error!("Optimize is expected to return one partition result");
  };
  let Some(sequences) = result.sequences.take() else {
    return make_internal_error!("Optimize returned no reconstructed sequences");
  };
  let maps = optimize_output_maps(&graph, sequences, UnknownMutationFilter::hiding_unknown(unknown))?;

  let topology_order = args
    .topology_order
    .resolve_topology_order(&graph, &names, None, args.tree_dialect.dialect())?;
  topology_order.apply(&mut graph, &names, &branch_lengths)?;
  stages.report("Writing output", 0.9, "");

  write_gtr_output(&resolved, models, "no GTR model was fitted", log)?;

  let trees = OptimizeTrees {
    graph: &graph,
    alphabet: &result.alphabet,
    names: &names,
    branch_lengths: &branch_lengths,
    maps: &maps,
    mutation_units: matches!(args.divergence_units, DivergenceUnits::Mutations),
    alignment: args.alignment.alignment.first().map(PathBuf::as_path),
    input_tree: args.tree(),
  };
  write_optimize_trees(&trees, &resolved, log)?;

  stages.report("Done", 1.0, "");

  Ok(())
}

fn write_optimize_trees(
  trees: &OptimizeTrees<'_>,
  resolved: &ResolvedOutputs,
  log: &dyn LogSink,
) -> Result<(), Report> {
  let annotated = AnnotatedGraph {
    graph: trees.graph,
    names: trees.names,
    divergence_branch_lengths: trees.branch_lengths,
    time_branch_lengths: None,
    divergence: Divergence::CumulativeBranchLength,
    sequences: Some(TreeSequences {
      alphabet: trees.alphabet,
      root_sequence: &trees.maps.root_sequence,
      edge_mutations: &trees.maps.edge_mutations,
      mutation_counts: trees.mutation_units.then_some(&trees.maps.edge_mutation_counts),
      amino_acids: None,
    }),
    dates: None,
    traits: None,
  };
  write_graph_outputs(&annotated, &resolved.tree_outputs)?;
  let Some(tree) = tree_view_for_outputs(&annotated, resolved)? else {
    return Ok(());
  };
  write_tree_outputs(&tree, &resolved.tree_outputs, CommandKind::Optimize, log)?;
  if let Some(path) = resolved.path(OutputSelection::AugurNodeData) {
    let run = RefineRun {
      alignment: trees.alignment,
      input_tree: Some(trees.input_tree),
      clock_model: None,
    };
    write_augur_node_data_refine(&tree, &run, path)?;
    progress_info!(log, "Wrote augur node data JSON to {path}", path = path.display());
  }
  Ok(())
}

struct OptimizeTrees<'a> {
  graph: &'a Graph,
  alphabet: &'a Alphabet,
  names: &'a BTreeMap<GraphNodeKey, Option<String>>,
  branch_lengths: &'a BTreeMap<GraphEdgeKey, Option<f64>>,
  maps: &'a OptimizeOutputMaps,
  mutation_units: bool,
  alignment: Option<&'a Path>,
  input_tree: &'a Path,
}

fn optimize_output_maps(
  graph: &Graph,
  sequences: PartitionSequences,
  filter: UnknownMutationFilter,
) -> Result<OptimizeOutputMaps, Report> {
  let PartitionSequences {
    root_sequence,
    edge_mutations,
    edge_mutation_counts,
  } = sequences;
  Ok(OptimizeOutputMaps {
    root_sequence,
    edge_mutations: filter.reported_edge_mutations(graph, edge_mutations)?,
    edge_mutation_counts,
  })
}

struct OptimizeOutputMaps {
  root_sequence: Seq,
  edge_mutations: BTreeMap<GraphEdgeKey, Vec<Mutation>>,
  edge_mutation_counts: BTreeMap<GraphEdgeKey, usize>,
}
