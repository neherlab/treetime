use crate::alphabet::alphabet::Alphabet;
use crate::commands::ancestral::fitch::{compress_sequences, get_common_length};
use crate::commands::ancestral::marginal::{initialize_marginal, update_marginal};
use crate::commands::optimize::args::TreetimeOptimizeArgs;
use crate::commands::optimize::optimize_unified::{initial_guess_mixed, run_optimize_mixed};
use crate::commands::optimize::partition_ops::{PartitionOptimizeOps, PartitionOptimizeVec};
use crate::gtr::get_gtr::{GtrModelName, JC69Params, jc69, log_gtr, write_gtr_json};
use crate::representation::algo::infer_dense::infer_dense;
use crate::representation::partition::marginal_dense::PartitionMarginalDense;
use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
use crate::representation::payload::ancestral::GraphAncestral;
use eyre::Report;
use itertools::Itertools;
use itertools::chain;
use log::debug;
use maplit::btreemap;
use parking_lot::RwLock;
use serde::Serialize;
use std::path::Path;
use std::sync::Arc;
use treetime_graph::edge::GraphEdge;
use treetime_graph::node::GraphNode;
use treetime_io::fasta::read_many_fasta;
use treetime_io::nex::{NexWriteOptions, nex_write_file};
use treetime_io::nwk::{EdgeToNwk, NodeToNwk, NwkWriteOptions, nwk_read_file, nwk_write_file};
use treetime_utils::fmt::float::float_to_significant_digits;

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
  } = args;

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
    // For now, use JC69 for both partitions to avoid GTR inference issues with dense partitions
    for partition in &sparse_partitions {
      let gtr = jc69(JC69Params::default())?; // FIXME: allow other models and model inference
      log_gtr(&gtr, GtrModelName::JC69);
      write_gtr_json(&gtr, GtrModelName::JC69, outdir)?;
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

    for partition in &dense_partitions {
      let gtr = jc69(JC69Params::default())?; // FIXME: allow other models and model inference
      log_gtr(&gtr, GtrModelName::JC69);
      partition.write_arc().gtr = gtr;
    }

    dense_partitions
  };

  // Run marginal reconstruction once to initialize edge data before optimization
  update_marginal(&graph, &sparse_partitions)?;
  initialize_marginal(&graph, &dense_partitions, &aln)?;

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

    run_optimize_mixed(&graph, &mixed_partitions)?;
    lh_prev = total_lh;
  }

  write_graph(outdir, &graph)?;
  Ok(())
}

fn collect_optimize_partitions(
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

fn write_graph<N, E, D>(outdir: impl AsRef<Path>, graph: &treetime_graph::graph::Graph<N, E, D>) -> Result<(), Report>
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
