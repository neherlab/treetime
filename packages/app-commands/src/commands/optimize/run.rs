use crate::commands::optimize::args::TreetimeOptimizeArgs;
use crate::commands::shared::gtr_output::{partition_models, write_gtr_output};
use crate::commands::shared::partition::log_partitions;
use crate::commands::shared::sequence_inputs::SequenceSources;
use crate::commands::shared::tree_input::read_input_tree;
use app_output::annotated_graph::{AnnotatedGraph, Divergence, TreeSequences};
use app_output::augur_node_data_refine::{RefineRun, total_mutation_counts, write_augur_node_data_refine};
use app_output::command_outputs::CommandKind;
use app_output::mutation_filter::UnknownMutationFilter;
use app_output::output_plan::{OutputSelection, ResolvedOutputs};
use app_output::tree_output::{write_run_tree_outputs, write_tree_outputs_per_alignment};
use eyre::Report;
use std::collections::BTreeMap;
use treetime::cancel::Cancel;
use treetime::optimize::pipeline::{self, OptimizeInput, OptimizeParams};
use treetime::partition::result::AlignmentSequences;
use treetime::partition::split::SequenceInput;
use treetime::progress::{LogSink, StageSink};
use treetime::progress_info;
use treetime::seq::mutation::Mutation;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::EdgeMap;
use treetime_graph::node_names::NodeNames;
use treetime_primitives::Seq;
use treetime_utils::make_internal_error;

pub fn run_optimize(
  args: &TreetimeOptimizeArgs,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<(), Report> {
  let resolved = args.outputs.resolve()?;
  cancel.check()?;
  stages.report("Reading input", 0.0, "");

  let nwk_parsed = read_input_tree(&args.tree, args.newick_dialect, log)?;
  let names = nwk_parsed.names;
  let graph = nwk_parsed.graph;
  let branch_lengths = nwk_parsed.branch_lengths;
  let SequenceSources { alignments, partitions } = args.sequences.read_partitions(&args.model, &graph, &names, log)?;
  log_partitions(&partitions, log);
  let sequences = SequenceInput {
    alignments: alignments.into_iter().map(|source| source.source).collect(),
    partitions,
  };

  let params = OptimizeParams {
    dense: args.dense,
    initial_guess: args.initial_branch_lengths,
    reroot_spec: args
      .reroot
      .as_ref()
      .map(|reroot| reroot.spec(&graph, &names))
      .transpose()?,
    loop_settings: args.loop_settings,
  };

  let input = OptimizeInput {
    graph,
    names,
    sequences,
    branch_lengths,
  };

  let output = pipeline::run(&params, input, cancel, stages, log).map_err(|err| err.into_report())?;
  let pipeline::OptimizeOutput {
    mut graph,
    alignments: mut results,
    partitions,
    branch_lengths,
    names,
  } = output;
  let models = partition_models(&partitions, &results, None)?;
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

  let topology_order = args.topology_order.resolve(&graph, &names, None, args.newick_dialect)?;
  topology_order.apply(&mut graph, &names, &branch_lengths)?;
  stages.report("Writing output", 0.9, "");

  write_gtr_output(&resolved, &models, "no GTR model was fitted", log)?;

  let mutation_counts = args
    .mutation_units
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
        comment_positions: args.positions,
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
      alignments: &args.sequences.alignments,
      mutation_counts: mutation_counts.as_ref(),
      input_tree: Some(args.tree.as_path()),
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
  write_tree_outputs_per_alignment(&annotated, trees.sequences, resolved, CommandKind::Optimize, log)
}

struct OptimizeTrees<'a> {
  graph: &'a Graph,
  names: &'a NodeNames,
  branch_lengths: &'a EdgeMap<Option<f64>>,
  sequences: &'a BTreeMap<usize, TreeSequences<'a>>,
  run: RefineRun<'a>,
}

fn optimize_output_maps(
  graph: &Graph,
  sequences: AlignmentSequences,
  filter: UnknownMutationFilter,
) -> Result<OptimizeOutputMaps, Report> {
  let AlignmentSequences {
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
  edge_mutations: EdgeMap<Vec<Mutation>>,
  edge_mutation_counts: EdgeMap<usize>,
}
