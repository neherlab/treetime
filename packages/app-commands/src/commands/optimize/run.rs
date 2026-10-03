use crate::commands::optimize::args::TreetimeOptimizeArgs;
use crate::commands::shared::output_args::DivergenceUnits;
use crate::commands::shared::resolve_outputs::ResolveOutputs;
use app_output::EdgeMutationCommentProvider;
use app_output::augur_node_data_optimize::write_augur_node_data_json;
use app_output::gtr::write_gtr_json;
use app_output::mutation_filter::UnknownMutationFilter;
use app_output::optimize_result::{OptimizeNodeOut, OptimizeOutputMaps};
use app_output::optimize_tree_output::write_optimize_tree_outputs;
use app_output::output_plan::OutputSelection;
use eyre::Report;
use std::collections::BTreeMap;
use std::path::PathBuf;
use treetime::alphabet::alphabet::Alphabet;
use treetime::cancel::Cancel;
use treetime::gtr::get_gtr::GtrOutput;
use treetime::optimize::pipeline::{self, OptimizeInput, OptimizeParams};
use treetime::partition::marginal::reconstruction::MarginalReconstruction;
use treetime::progress::{LogSink, StageSink};
use treetime::progress_info;
use treetime::seq::gap_fill::apply_gap_fill;
use treetime::seq::mutation::{MutationTrack, edge_state_change_counts};
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::fasta::read_many_fasta_path;
use treetime_io::nwk::CommentProviders;
use treetime_io::nwk::nwk_read_file;
use treetime_primitives::AlignmentRecord;

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
  let mut aln = read_many_fasta_path(&args.alignment.alignment, &alphabet)?;
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

  let nodes: BTreeMap<GraphNodeKey, OptimizeNodeOut> = graph
    .get_nodes()
    .map(|node| {
      let key = node.key();
      (
        key,
        OptimizeNodeOut {
          name: names[&key].clone(),
          branch_support: confidences.get(&key).copied().flatten(),
        },
      )
    })
    .collect();

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::Gtr) {
    let gtr_output = GtrOutput::builder().gtr(&gtr).model_name(model_name).build();
    write_gtr_json(&gtr_output, path)?;
  }

  if !resolved.tree_outputs.is_empty() {
    let mutation_provider = EdgeMutationCommentProvider::new(&maps.edge_mutations, &graph);
    let providers = CommentProviders::new().with(&mutation_provider);
    write_optimize_tree_outputs(
      &graph,
      &nodes,
      &branch_lengths,
      &maps,
      &resolved.tree_outputs,
      &providers,
    )?;
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::AugurNodeData) {
    let mutation_counts = match args.divergence_units {
      DivergenceUnits::Mutations => Some(&maps.edge_mutation_counts),
      DivergenceUnits::MutationsPerSite => None,
    };

    let alignment = args.alignment.alignment.first().map(PathBuf::as_path);
    write_augur_node_data_json(
      &graph,
      &nodes,
      &branch_lengths,
      alignment,
      Some(args.tree()),
      mutation_counts,
      path,
    )?;
    progress_info!(log, "Wrote augur node data JSON to {path}", path = path.display());
  }

  stages.report("Done", 1.0, "");

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
