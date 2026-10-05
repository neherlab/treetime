use crate::commands::optimize::args::TreetimeOptimizeArgs;
use crate::commands::shared::alignment::read_alignment;
use crate::commands::shared::output_args::DivergenceUnits;
use crate::commands::shared::resolve_outputs::ResolveOutputs;
use app_output::annotated_graph::{AnnotatedGraph, Divergence, TreeSequences};
use app_output::augur_node_data_refine::{RefineRun, write_augur_node_data_refine};
use app_output::mutation_filter::UnknownMutationFilter;
use app_output::output_plan::{CommandKind, OutputSelection, ResolvedOutputs};
use app_output::tree_output::{tree_view_for_outputs, write_graph_outputs, write_tree_outputs};
use eyre::Report;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use treetime::alphabet::alphabet::Alphabet;
use treetime::cancel::Cancel;
use treetime::gtr::get_gtr::GtrOutput;
use treetime::optimize::pipeline::{self, OptimizeInput, OptimizeParams};
use treetime::partition::marginal::reconstruction::MarginalReconstruction;
use treetime::progress::{LogSink, StageSink};
use treetime::progress_info;
use treetime::seq::gap_fill::apply_gap_fill;
use treetime::seq::mutation::{Mutation, MutationTrack, edge_state_change_counts};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::nwk::nwk_read_file;
use treetime_primitives::{AlignmentRecord, Seq};
use treetime_utils::io::json::{JsonPretty, json_write_file};

pub fn run_optimize(
  args: &TreetimeOptimizeArgs,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<(), Report> {
  cancel.check()?;
  stages.report("Reading input", 0.0, "");

  let alphabet = Alphabet::new(args.alphabet_args.alphabet_name().unwrap_or_default())?;
  let gap_fill = args.gap_fill_args.effective_gap_fill();
  let mut aln = read_alignment(&args.alignment.alignment, &alphabet)?;
  for record in &mut aln {
    apply_gap_fill(&mut record.seq, gap_fill, alphabet.gap(), alphabet.unknown());
  }
  let nwk_parsed = nwk_read_file(args.tree())?;
  let confidences = nwk_parsed.confidences();
  let names = nwk_parsed.names();
  let graph = nwk_parsed.graph;
  let branch_lengths = nwk_parsed.branch_lengths;

  let resolved = args.resolve_outputs()?;

  let params = OptimizeParams {
    model: args.model_args.model_name(),
    dense: args.dense,
    max_iter: args.max_iter,
    dp: args.dp,
    damping: args.damping,
    opt_method: args.opt_method,
    initial_guess: args.branch_length_initial_guess,
    no_indels: args.no_indels,
    reroot_spec: args.reroot_spec(),
    topology_ops: args.topology_ops,
  };

  let unknown = alphabet.unknown();
  let input = OptimizeInput {
    graph,
    alphabet,
    sequences: aln.into_iter().map(AlignmentRecord::from).collect(),
    branch_lengths,
  };

  let output = pipeline::run(&params, input, &names, cancel, stages, log).map_err(|err| err.into_report())?;
  let pipeline::OptimizeOutput {
    mut graph,
    gtr,
    model_name,
    reconstruction,
    branch_lengths,
    names,
  } = output;

  let maps = gather_optimize_output_maps(&graph, &reconstruction, UnknownMutationFilter::hiding_unknown(unknown))?;

  let topology_order = args.topology_order.resolve_topology_order(&graph, &names, None)?;
  topology_order.apply(&mut graph, &names, &branch_lengths)?;
  stages.report("Writing output", 0.9, "");

  if let Some(path) = resolved.path(OutputSelection::Gtr) {
    let gtr_output = GtrOutput::builder().gtr(&gtr).model_name(model_name).build();
    json_write_file(path, &gtr_output, JsonPretty(true))?;
  }

  let trees = OptimizeTrees {
    graph: &graph,
    names: &names,
    branch_lengths: &branch_lengths,
    branch_support: &confidences,
    maps: &maps,
    mutation_units: matches!(args.divergence_units, DivergenceUnits::Mutations),
    alignment: args.alignment.alignment.first().map(PathBuf::as_path),
    input_tree: args.tree(),
  };
  write_optimize_trees(&trees, &resolved, log)?;

  stages.report("Done", 1.0, "");

  Ok(())
}

struct OptimizeTrees<'a> {
  graph: &'a Graph,
  names: &'a BTreeMap<GraphNodeKey, Option<String>>,
  branch_lengths: &'a BTreeMap<GraphEdgeKey, Option<f64>>,
  branch_support: &'a BTreeMap<GraphNodeKey, Option<f64>>,
  maps: &'a OptimizeOutputMaps,
  mutation_units: bool,
  alignment: Option<&'a Path>,
  input_tree: &'a Path,
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
    branch_support: Some(trees.branch_support),
    sequences: Some(TreeSequences {
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
      branch_support: Some(trees.branch_support),
    };
    write_augur_node_data_refine(&tree, &run, path)?;
    progress_info!(log, "Wrote augur node data JSON to {path}", path = path.display());
  }
  Ok(())
}

fn gather_optimize_output_maps(
  graph: &Graph,
  reconstruction: &MarginalReconstruction,
  filter: UnknownMutationFilter,
) -> Result<OptimizeOutputMaps, Report> {
  let edge_mutations = graph
    .get_edges()
    .map(|edge| {
      let key = edge.key();
      Ok((
        key,
        reconstruction.edge_mutations(graph, key, &MutationTrack::Nucleotide)?,
      ))
    })
    .collect::<Result<BTreeMap<_, _>, Report>>()?;
  let edge_mutation_counts = edge_state_change_counts(&edge_mutations, reconstruction.alphabet())?;
  Ok(OptimizeOutputMaps {
    root_sequence: reconstruction.root_sequence(graph)?,
    edge_mutations: filter.reported_edge_mutations(graph, edge_mutations)?,
    edge_mutation_counts,
  })
}

struct OptimizeOutputMaps {
  root_sequence: Seq,
  edge_mutations: BTreeMap<GraphEdgeKey, Vec<Mutation>>,
  edge_mutation_counts: BTreeMap<GraphEdgeKey, usize>,
}
