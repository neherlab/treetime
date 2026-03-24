use crate::alphabet::alphabet::Alphabet;
use crate::commands::ancestral::fitch::{compress_sequences, get_common_length};
use crate::commands::ancestral::marginal::{initialize_marginal, update_marginal};
use crate::commands::optimize::args::TreetimeOptimizeArgs;
use crate::commands::optimize::optimize_unified::{initial_guess_mixed, run_optimize_mixed};
use crate::commands::optimize::partition_ops::{PartitionOptimizeOps, PartitionOptimizeVec};
use crate::gtr::get_gtr::{JC69Params, get_gtr_dense, get_gtr_sparse, jc69, write_gtr_json};
use crate::representation::algo::infer_dense::infer_dense;
use crate::representation::partition::marginal_dense::PartitionMarginalDense;
use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
use crate::representation::payload::ancestral::GraphAncestral;
use eyre::Report;
use itertools::{Itertools, chain, izip};
use log::debug;
use maplit::btreemap;
use num_traits::pow::pow;
use parking_lot::RwLock;
use serde::Serialize;
use std::path::Path;
use std::sync::Arc;
use treetime_graph::edge::{GraphEdge, HasBranchLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNode;
use treetime_io::fasta::read_many_fasta;
use treetime_io::nex::{NexWriteOptions, nex_write_file};
use treetime_io::nwk::{EdgeToNwk, NodeToNwk, NwkWriteOptions, nwk_read_file, nwk_write_file};
use treetime_utils::fmt::float::float_to_significant_digits;
use treetime_utils::make_error;

#[derive(Clone, Debug, Default)]
pub struct TreetimeOptimizeParams {
  pub sample_from_profile: bool,
  pub fixed_pi: bool,
}

pub fn run_optimize(args: &TreetimeOptimizeArgs) -> Result<(), Report> {
  let TreetimeOptimizeArgs {
    input_fastas,
    tree,
    alphabet,
    model_name,
    dense,
    outdir,
    max_iter,
    dp,
    damping,
  } = args;

  if !(0.0..1.0).contains(damping) {
    return make_error!("--damping must be in [0.0, 1.0), got {damping}");
  }

  // FIXME: For demonstration of mixed partitions, we create a sparse and a dense partition,
  // both containing full sequence. Make it configurable.

  // TODO: currently unused, but we could implement automatic suggestion of dense/sparse mode for individual partitions
  let _dense = dense.unwrap_or_else(infer_dense);

  let alphabet_sparse = Alphabet::new(alphabet.unwrap_or_default())?;
  let alphabet_dense = Alphabet::new(alphabet.unwrap_or_default())?;

  let aln = read_many_fasta(input_fastas, &alphabet_sparse)?;

  let graph: GraphAncestral = nwk_read_file(tree)?;

  let sparse_partitions = {
    #[allow(clippy::iter_on_single_items)]
    let sparse_partitions = [PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?, // FIXME: dummy temporary gtr should not be needed here
      alphabet: alphabet_sparse,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }]
    .into_iter()
    .map(|p| Arc::new(RwLock::new(p)))
    .collect_vec();

    compress_sequences(&graph, &sparse_partitions, &aln)?;

    // FIXME: chicken & egg problem: to get a gtr we need partitions, to get partitions we need a gtr
    // FIXME: spaghetti code: dummy gtr is replaced by real gtr here
    for partition in &sparse_partitions {
      let gtr = get_gtr_sparse(model_name, partition, &graph)?;
      write_gtr_json(&gtr, *model_name, outdir, Some("sparse"))?;
      partition.write_arc().gtr = gtr;
    }

    sparse_partitions
  };

  let dense_partitions = {
    #[allow(clippy::iter_on_single_items)]
    let dense_partitions = [PartitionMarginalDense {
      index: 1,
      gtr: jc69(JC69Params::default())?, // FIXME: dummy temporary gtr should not be needed here
      alphabet: alphabet_dense,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }]
    .into_iter()
    .map(|p| Arc::new(RwLock::new(p)))
    .collect_vec();

    dense_partitions
  };

  // Run marginal reconstruction to initialize edge data before optimization
  update_marginal(&graph, &sparse_partitions)?;
  initialize_marginal(&graph, &dense_partitions, &aln)?;

  // FIXME: chicken & egg problem: to get a gtr we need partitions, to get partitions we need a gtr
  // FIXME: spaghetti code: dummy gtr is replaced by real gtr here
  // Dense GTR requires profiles populated via initialize_marginal + update_marginal.
  // Run first marginal pass with dummy GTR, then replace with real model.
  update_marginal(&graph, &dense_partitions)?;
  for partition in &dense_partitions {
    let gtr = get_gtr_dense(model_name, partition, &graph)?;
    write_gtr_json(&gtr, *model_name, outdir, Some("dense"))?;
    partition.write_arc().gtr = gtr;
  }

  // Re-run marginal with inferred GTR so node posteriors used by
  // initial_guess_mixed reflect the real model, not the dummy JC69.
  update_marginal(&graph, &dense_partitions)?;

  let mixed_partitions = collect_optimize_partitions(&dense_partitions, &sparse_partitions);

  initial_guess_mixed(&graph, &mixed_partitions)?;

  let mut lh_prev = f64::MIN;
  for i in 0..*max_iter {
    let sparse_lh = update_marginal(&graph, &sparse_partitions)?;
    let dense_lh = update_marginal(&graph, &dense_partitions)?;
    let total_lh = sparse_lh + dense_lh;

    debug!(
      "Iteration {}: likelihood {} (sparse: {}, dense: {})",
      i + 1,
      float_to_significant_digits(total_lh, 7),
      float_to_significant_digits(sparse_lh, 7),
      float_to_significant_digits(dense_lh, 7)
    );

    if (total_lh - lh_prev).abs() < dp.abs() {
      break;
    }

    let old_branch_lengths = save_branch_lengths(&graph);
    run_optimize_mixed(&graph, &mixed_partitions)?;
    apply_damping(&graph, &old_branch_lengths, *damping, i);
    lh_prev = total_lh;
  }

  write_graph(outdir, &graph)?;
  Ok(())
}

pub(crate) fn collect_optimize_partitions(
  dense_partitions: &[Arc<RwLock<PartitionMarginalDense>>],
  sparse_partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
) -> PartitionOptimizeVec {
  chain!(
    dense_partitions
      .iter()
      .cloned()
      .map(|partition| -> Arc<RwLock<dyn PartitionOptimizeOps>> { partition }),
    sparse_partitions
      .iter()
      .cloned()
      .map(|partition| -> Arc<RwLock<dyn PartitionOptimizeOps>> { partition })
  )
  .collect_vec()
}

/// Save current branch lengths for all edges in graph traversal order.
pub(crate) fn save_branch_lengths<N, E, D>(graph: &Graph<N, E, D>) -> Vec<f64>
where
  N: GraphNode,
  E: GraphEdge + HasBranchLength,
  D: Send + Sync,
{
  graph
    .get_edges()
    .iter()
    .map(|edge_ref| edge_ref.read_arc().payload().read_arc().branch_length().unwrap_or(0.0))
    .collect_vec()
}

/// Blend optimized branch lengths with saved old values using exponential damping.
///
/// At iteration `i` (0-based), each branch length becomes:
///   bl = bl_optimized * (1 - damping^(i+1)) + bl_old * damping^(i+1)
///
/// When `damping == 0.0`, damping_factor = 0 and the optimized value is kept unchanged.
/// Early iterations take conservative steps; later iterations approach the full update.
pub(crate) fn apply_damping<N, E, D>(graph: &Graph<N, E, D>, old_branch_lengths: &[f64], damping: f64, iteration: usize)
where
  N: GraphNode,
  E: GraphEdge + HasBranchLength,
  D: Send + Sync,
{
  if damping == 0.0 {
    return;
  }
  debug_assert_eq!(
    graph.get_edges().len(),
    old_branch_lengths.len(),
    "apply_damping: edge count changed between save and apply"
  );
  let damping_factor = pow(damping, iteration + 1);
  let new_weight = 1.0 - damping_factor;
  for (edge_ref, &old_bl) in izip!(graph.get_edges(), old_branch_lengths) {
    let mut edge = edge_ref.write_arc().payload().write_arc();
    let optimized_bl = edge.branch_length().unwrap_or(0.0);
    let damped_bl = optimized_bl * new_weight + old_bl * damping_factor;
    edge.set_branch_length(Some(damped_bl));
  }
}

fn write_graph<N, E, D>(outdir: impl AsRef<Path>, graph: &Graph<N, E, D>) -> Result<(), Report>
where
  N: GraphNode + NodeToNwk + Serialize,
  E: GraphEdge + EdgeToNwk + Serialize,
  D: Send + Sync + Default + Serialize,
{
  // json_write_file(
  //   outdir.as_ref().join("annotated_tree.graph.json"),
  //   &graph,
  //   JsonPretty(true),
  // )?;

  nwk_write_file(
    outdir.as_ref().join("annotated_tree.nwk"),
    graph,
    &NwkWriteOptions::default(),
  )?;

  nex_write_file(
    outdir.as_ref().join("annotated_tree.nexus"),
    graph,
    &NexWriteOptions::default(),
  )?;

  Ok(())
}
