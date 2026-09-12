use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::pipeline::{DenseReconstruction, SparseReconstruction};
use crate::commands::optimize::args::TreetimeOptimizeArgs;
use crate::commands::optimize::augur_node_data::write_augur_node_data_json;
use crate::commands::optimize::result::{EdgeOut, OptimizeNodeOut, OptimizeOutputMaps, OptimizeResult};
use crate::commands::shared::output::{DivergenceUnits, OutputSelection};
use crate::commands::shared::resolve_outputs::ResolveOutputs;
use crate::commands::shared::tree_output::write_optimize_tree_outputs;
use crate::gtr::get_gtr::{GtrOutput, write_gtr_json};
use crate::make_error;
use crate::optimize::pipeline::{self, OptimizeInput, OptimizeParams};
use crate::partition::traits::{EdgeMutationCommentProvider, PartitionBranchOps};
use crate::seq::gap_fill::apply_gap_fill;
use crate::seq::mutation::MutationTrack;
use eyre::Report;
use log::info;
use std::collections::BTreeMap;
use std::path::PathBuf;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::fasta::read_many_fasta;
use treetime_io::nwk::CommentProviders;
use treetime_io::nwk::{NwkParse, nwk_read_file};

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
  let NwkParse {
    graph,
    confidences,
    names,
    branch_lengths,
  } = nwk_read_file(args.tree())?;

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
    branch_lengths,
  };

  let output = pipeline::run(&params, input, &names, progress)?;
  let pipeline::OptimizeOutput {
    mut graph,
    gtr,
    model_name,
    sparse_partitions,
    dense_partitions,
    branch_lengths,
    names,
  } = output;

  // Gather the per-node/per-edge sequence and mutation values off the pipeline-local partitions into
  // plain value maps the output writers consume. This is the only place that reads sequences and
  // mutations from the partition; the tree, node-data, and Newick-comment writers read the maps
  // instead. Node and edge keys stay stable through topology ordering, so gathering before it is
  // bit-identical.
  let maps = gather_optimize_output_maps(&graph, &sparse_partitions, &dense_partitions)?;
  let has_partitions = !sparse_partitions.is_empty() || !dense_partitions.is_empty();

  let topology_order = args.topology_order.resolve_topology_order(&graph, &names, None)?;
  topology_order.apply(&mut graph, &names, &branch_lengths)?;
  progress.report("Writing output", 0.9, "");

  // Gather the per-node name/confidence into a keyed value map the output writers consume. The name
  // comes from the pipeline's post-loop name map (`names`); topology ordering only permutes
  // children, so the map still matches the ordered tree. The optimized per-edge branch lengths come
  // from the loop result (`branch_lengths`); the writers read sequences from the gathered maps.
  let nodes: BTreeMap<GraphNodeKey, OptimizeNodeOut> = graph
    .get_nodes()
    .iter()
    .map(|node| {
      let node = node.read_arc();
      let key = node.key();
      (
        key,
        OptimizeNodeOut {
          name: names[&key].clone(),
          confidence: confidences.get(&key).copied().flatten(),
        },
      )
    })
    .collect();
  let edges: BTreeMap<GraphEdgeKey, EdgeOut> = branch_lengths
    .iter()
    .map(|(&key, &branch_length)| (key, EdgeOut { branch_length }))
    .collect();

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::Gtr) {
    let gtr_output = GtrOutput::new(&gtr, model_name);
    write_gtr_json(&gtr_output, path)?;
  }

  if !resolved.tree_outputs.is_empty() {
    // Dense and sparse reconstructions annotate Newick/Nexus nodes with their inbound mutations; the
    // partition-less case emits no such comments. The comment provider now reads the gathered per-edge
    // mutation map rather than the partition.
    if maps.root_sequence.is_some() {
      let provider = EdgeMutationCommentProvider::new(&maps.edge_mutations, &graph);
      let providers = CommentProviders::new().with(&provider);
      write_optimize_tree_outputs(
        &graph,
        &nodes,
        &branch_lengths,
        &maps,
        &resolved.tree_outputs,
        &providers,
      )?;
    } else {
      write_optimize_tree_outputs(
        &graph,
        &nodes,
        &branch_lengths,
        &maps,
        &resolved.tree_outputs,
        &CommentProviders::new(),
      )?;
    }
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::AugurNodeData) {
    let mutation_counts = match args.divergence_units {
      DivergenceUnits::Mutations => {
        if !has_partitions {
          return make_error!(
            "--divergence-units=mutations requires ancestral reconstruction but no partitions are available"
          );
        }
        Some(
          maps
            .edge_subs
            .iter()
            .map(|(&key, subs)| (key, subs.len()))
            .collect::<BTreeMap<GraphEdgeKey, usize>>(),
        )
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

  Ok(OptimizeResult { graph, nodes, edges })
}

/// Gather the per-node nucleotide sequences, root sequence, per-edge nucleotide mutations, and per-edge
/// substitutions the output writers read off the optimize partition.
pub(crate) fn gather_optimize_output_maps(
  graph: &Graph,
  sparse_partitions: &[SparseReconstruction],
  dense_partitions: &[DenseReconstruction],
) -> Result<OptimizeOutputMaps, Report> {
  if let Some(family) = dense_partitions.first() {
    gather_optimize_partition_maps(graph, &family.readout())
  } else if let Some(family) = sparse_partitions.first() {
    gather_optimize_partition_maps(graph, &family.readout())
  } else {
    Ok(OptimizeOutputMaps::default())
  }
}

fn gather_optimize_partition_maps(
  graph: &Graph,
  partition: &dyn PartitionBranchOps,
) -> Result<OptimizeOutputMaps, Report> {
  let root_sequence = Some(partition.root_sequence(graph)?);
  let node_sequences = graph
    .get_nodes()
    .iter()
    .map(|node| {
      let key = node.read_arc().key();
      (key, partition.node_sequence(key))
    })
    .collect();
  let edge_mutations = graph
    .get_edges()
    .iter()
    .map(|edge| {
      let key = edge.read_arc().key();
      Ok((key, partition.edge_mutations(graph, key, MutationTrack::Nucleotide)?))
    })
    .collect::<Result<BTreeMap<_, _>, Report>>()?;
  let edge_subs = graph
    .get_edges()
    .iter()
    .map(|edge| {
      let key = edge.read_arc().key();
      Ok((key, partition.edge_subs(graph, key)?))
    })
    .collect::<Result<BTreeMap<_, _>, Report>>()?;
  Ok(OptimizeOutputMaps {
    root_sequence,
    node_sequences,
    edge_mutations,
    edge_subs,
  })
}
