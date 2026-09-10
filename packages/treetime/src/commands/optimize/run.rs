use crate::alphabet::alphabet::Alphabet;
use crate::commands::optimize::args::TreetimeOptimizeArgs;
use crate::commands::optimize::augur_node_data::write_augur_node_data_json;
use crate::commands::optimize::result::{EdgeOut, OptimizeGraphData, OptimizeNodeOut, OptimizeResult};
use crate::commands::shared::output::{DivergenceUnits, OutputSelection};
use crate::commands::shared::resolve_outputs::ResolveOutputs;
use crate::commands::shared::tree_output::write_optimize_tree_outputs;
use crate::gtr::get_gtr::{GtrOutput, write_gtr_json};
use crate::make_error;
use crate::optimize::pipeline::{self, OptimizeInput, OptimizeParams};
use crate::partition::traits::MutationCommentProvider;
use crate::seq::div::compute_edge_mutation_counts;
use crate::seq::gap_fill::apply_gap_fill;
use eyre::Report;
use log::info;
use std::collections::BTreeMap;
use std::path::PathBuf;
use treetime_graph::edge::{GraphEdgeKey, HasBranchLength};
use treetime_graph::node::GraphNodeKey;
use treetime_io::fasta::read_many_fasta;
use treetime_io::nwk::CommentProviders;
use treetime_io::nwk::nwk_read_file;

pub fn run_optimize(
  args: &TreetimeOptimizeArgs,
  progress: &dyn crate::progress::ProgressSink,
) -> Result<OptimizeResult, Report> {
  progress.check_cancelled()?;
  progress.report("Reading input", 0.0, "");

  let alphabet = Alphabet::new(args.alphabet_args.alphabet.unwrap_or_default())?;
  let gap_fill = args.gap_fill_args.effective_gap_fill();
  let mut aln = read_many_fasta(&args.alignment.alignment, &alphabet)?;
  for record in &mut aln {
    apply_gap_fill(&mut record.seq, gap_fill, alphabet.gap(), alphabet.unknown());
  }
  let graph = nwk_read_file(args.tree())?;

  let resolved = args.resolve_outputs()?;

  let params = OptimizeParams {
    model: args.model_args.model,
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

  let input = OptimizeInput {
    graph,
    alphabet,
    sequences: aln,
  };

  let output = pipeline::run(&params, input, progress)?;
  let pipeline::OptimizeOutput {
    graph,
    gtr,
    model_name,
    sparse_partitions,
    dense_partitions,
  } = output;
  let mut graph = graph.map_data(OptimizeGraphData::new(
    gtr,
    model_name,
    sparse_partitions,
    dense_partitions,
  ));
  let topology_order = args.topology_order.resolve_topology_order(&graph, None)?;
  topology_order.apply(&mut graph)?;
  progress.report("Writing output", 0.9, "");

  // Gather the per-node name/confidence and the optimized per-edge branch length off the ordered
  // tree into keyed value maps the output writers consume. The writers still read sequences and
  // model metadata from the graph data slot; these maps carry the name, input-branch-support, and
  // branch-length reads that move off the payload.
  let nodes: BTreeMap<GraphNodeKey, OptimizeNodeOut> = graph
    .get_nodes()
    .iter()
    .map(|node| {
      let node = node.read_arc();
      let payload = node.payload().read_arc();
      (
        node.key(),
        OptimizeNodeOut {
          name: payload.name.clone(),
          confidence: payload.confidence,
        },
      )
    })
    .collect();
  let edges: BTreeMap<GraphEdgeKey, EdgeOut> = graph
    .get_edges()
    .iter()
    .map(|edge_ref| {
      let edge = edge_ref.read_arc();
      let branch_length = edge.payload().read_arc().branch_length();
      (edge.key(), EdgeOut { branch_length })
    })
    .collect();
  let branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>> =
    edges.iter().map(|(key, edge)| (*key, edge.branch_length)).collect();

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::Gtr) {
    let gtr_output = GtrOutput::new(&graph.data().gtr, graph.data().model_name);
    write_gtr_json(&gtr_output, path)?;
  }

  if !resolved.tree_outputs.is_empty() {
    if !graph.data().dense_partitions.is_empty() {
      let guard = graph.data().dense_partitions[0].read_arc();
      let provider = MutationCommentProvider::new(&*guard, &graph);
      let providers = CommentProviders::new().with(&provider);
      write_optimize_tree_outputs(&graph, &nodes, &branch_lengths, &resolved.tree_outputs, &providers)?;
    } else if !graph.data().sparse_partitions.is_empty() {
      let guard = graph.data().sparse_partitions[0].read_arc();
      let provider = MutationCommentProvider::new(&*guard, &graph);
      let providers = CommentProviders::new().with(&provider);
      write_optimize_tree_outputs(&graph, &nodes, &branch_lengths, &resolved.tree_outputs, &providers)?;
    } else {
      write_optimize_tree_outputs(&graph, &nodes, &branch_lengths, &resolved.tree_outputs, &CommentProviders::new())?;
    }
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::AugurNodeData) {
    let mutation_counts = match args.divergence_units {
      DivergenceUnits::Mutations => {
        let partition: &dyn crate::partition::traits::PartitionBranchOps = if !graph.data().dense_partitions.is_empty()
        {
          &*graph.data().dense_partitions[0].read_arc()
        } else if !graph.data().sparse_partitions.is_empty() {
          &*graph.data().sparse_partitions[0].read_arc()
        } else {
          return make_error!(
            "--divergence-units=mutations requires ancestral reconstruction but no partitions are available"
          );
        };
        Some(compute_edge_mutation_counts(&graph, partition)?)
      },
      DivergenceUnits::MutationsPerSite => None,
    };

    let alignment = args.alignment.alignment.first().map(PathBuf::as_path);
    write_augur_node_data_json(
      &graph,
      &nodes,
      &branch_lengths,
      alignment,
      Some(args.tree()),
      mutation_counts.as_ref(),
      path,
    )?;
    info!("Wrote augur node data JSON to {path}", path = path.display());
  }

  progress.report("Done", 1.0, "");

  let gtr = graph.data().gtr.clone();
  let model_name = graph.data().model_name;
  let sparse_partitions = graph.data().sparse_partitions.clone();
  let dense_partitions = graph.data().dense_partitions.clone();

  Ok(OptimizeResult {
    graph,
    nodes,
    edges,
    gtr,
    model_name,
    sparse_partitions,
    dense_partitions,
  })
}
