use crate::commands::optimize::args::TreetimeOptimizeArgs;
use crate::commands::shared::gtr_output::{partition_models, write_gtr_output};
use crate::commands::shared::output_args::DivergenceUnits;
use crate::commands::shared::output_args::ResolveOutputs;
use crate::commands::shared::sequence_inputs::{SourceSettings, read_sequence_sources};
use crate::commands::shared::tree_input::read_input_tree;
use app_output::annotated_graph::{AnnotatedGraph, Divergence, TreeSequences};
use app_output::augur_node_data_refine::{RefineRun, total_mutation_counts, write_augur_node_data_refine};
use app_output::mutation_filter::UnknownMutationFilter;
use app_output::output_plan::{CommandKind, OutputSelection, ResolvedOutputs};
use app_output::tree_output::{write_run_tree_outputs, write_tree_outputs_per_partition};
use eyre::Report;
use std::collections::BTreeMap;
use treetime::alphabet::alphabet::Alphabet;
use treetime::cancel::Cancel;
use treetime::optimize::params::OptimizeLoopSettings;
use treetime::optimize::pipeline::{self, OptimizeInput, OptimizeParams};
use treetime::partition::result::PartitionSequences;
use treetime::progress::{LogSink, StageSink};
use treetime::progress_info;
use treetime::seq::mutation::Mutation;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node_names::NodeNames;
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
  let nwk_parsed = read_input_tree(args.tree(), args.newick_dialect.dialect(), log)?;
  let names = nwk_parsed.names;
  let graph = nwk_parsed.graph;
  let branch_lengths = nwk_parsed.branch_lengths;
  let settings = SourceSettings {
    alphabet: &alphabet,
    gap_fill: Some(args.gap_fill_args.gap_fill()),
    missing_leaves: args.missing_sequences.limit(),
  };
  let sources = read_sequence_sources(&args.alignment, &settings, &graph, &names, log)?
    .into_iter()
    .map(|source| source.source)
    .collect();

  let resolved = args.resolve_outputs()?;

  let params = OptimizeParams {
    model: args.model.clone(),
    dense: args.representation.dense(),
    initial_guess: args.initial_branch_lengths,
    reroot_spec: args.reroot_spec(&graph, &names)?,
    loop_settings: OptimizeLoopSettings {
      max_iter: args.max_iter,
      dp: args.tolerance,
      damping: args.damping,
      opt_method: args.optimizer,
      indel_policy: args.indel_policy,
      topology_ops: args.topology_ops,
      partition_rates: args.partition_rates_args.partition_rates(),
    },
  };

  let input = OptimizeInput {
    graph,
    names,
    sources,
    branch_lengths,
  };

  let output = pipeline::run(&params, input, cancel, stages, log).map_err(|err| err.into_report())?;
  let pipeline::OptimizeOutput {
    mut graph,
    mut results,
    branch_lengths,
    names,
  } = output;
  let models = partition_models(&results, None);
  let maps = results
    .iter_mut()
    .map(|result| {
      let Some(sequences) = result.sequences.take() else {
        return make_internal_error!("Optimize returned no reconstructed sequences");
      };
      let filter = UnknownMutationFilter::hiding_unknown(result.alphabet.unknown());
      Ok((result.number, optimize_output_maps(&graph, sequences, filter)?))
    })
    .collect::<Result<BTreeMap<_, _>, Report>>()?;

  let topology_order =
    args
      .topology_order
      .resolve_topology_order(&graph, &names, None, args.newick_dialect.dialect())?;
  topology_order.apply(&mut graph, &names, &branch_lengths)?;
  stages.report("Writing output", 0.9, "");

  write_gtr_output(&resolved, models, "no GTR model was fitted", log)?;

  let mutation_counts = matches!(args.branch_lengths.divergence_units, DivergenceUnits::Mutations)
    .then(|| total_mutation_counts(maps.values().map(|maps| &maps.edge_mutation_counts)));
  let sequences = results
    .iter()
    .map(|result| {
      let maps = &maps[&result.number];
      let sequences = TreeSequences {
        alphabet: &result.alphabet,
        root_sequence: &maps.root_sequence,
        edge_mutations: &maps.edge_mutations,
        amino_acids: None,
        comment_positions: args.zero_based.positions(),
      };
      (result.number, sequences)
    })
    .collect();
  let trees = OptimizeTrees {
    graph: &graph,
    names: &names,
    branch_lengths: &branch_lengths,
    sequences: &sequences,
    run: RefineRun {
      alignments: &args.alignment.alignment,
      mutation_counts: mutation_counts.as_ref(),
      input_tree: Some(args.tree()),
      clock_model: None,
    },
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
    sequences: None,
    dates: None,
    traits: None,
  };
  let tree = write_run_tree_outputs(&annotated, resolved, CommandKind::Optimize, log)?;
  if let (Some(tree), Some(path)) = (&tree, resolved.path(OutputSelection::NodeData)) {
    write_augur_node_data_refine(tree, &trees.run, path)?;
    progress_info!(log, "Wrote augur node data JSON to {path}", path = path.display());
  }
  write_tree_outputs_per_partition(&annotated, trees.sequences, resolved, CommandKind::Optimize, log)
}

struct OptimizeTrees<'a> {
  graph: &'a Graph,
  names: &'a NodeNames,
  branch_lengths: &'a BTreeMap<GraphEdgeKey, Option<f64>>,
  sequences: &'a BTreeMap<usize, TreeSequences<'a>>,
  run: RefineRun<'a>,
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
