use crate::alphabet::alphabet::Alphabet;
use crate::commands::ancestral::fitch::{compress_sequences, get_common_length};
use crate::commands::ancestral::marginal_dense::run_marginal_dense;
use crate::commands::ancestral::marginal_sparse::run_marginal_sparse;
use crate::commands::optimize::args::TreetimeOptimizeArgs;
use crate::commands::optimize::optimize_dense::run_optimize_dense;
use crate::commands::optimize::optimize_sparse::{initial_guess_sparse, run_optimize_sparse};
use crate::graph::edge::GraphEdge;
use crate::graph::node::GraphNode;
use crate::gtr::get_gtr::{JC69Params, jc69};
use crate::io::fasta::read_many_fasta;
use crate::io::nex::{NexWriteOptions, nex_write_file};
use crate::io::nwk::{EdgeToNwk, NodeToNwk, NwkWriteOptions, nwk_read_file, nwk_write_file};
use crate::representation::graph_ancestral::GraphAncestral;
use crate::representation::infer_dense::infer_dense;
use crate::representation::partition_marginal_dense::PartitionMarginalDense;
use crate::representation::partition_marginal_sparse::PartitionMarginalSparse;
use crate::utils::float_fmt::float_to_significant_digits;
use eyre::Report;
use itertools::Itertools;
use log::debug;
use maplit::btreemap;
use parking_lot::RwLock;
use serde::Serialize;
use std::path::Path;
use std::sync::Arc;

// The initial guess for dense is not working well, but optimization works without revisit after settling on optimization algorithm
// use super::optimize_dense::initial_guess;

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

  let dense = dense.unwrap_or_else(infer_dense);

  let treat_gap_as_unknown = dense;
  let alphabet = Alphabet::new(alphabet.unwrap_or_default(), treat_gap_as_unknown)?;

  // TODO: avoid reading all sequences into memory somehow?
  let aln = read_many_fasta(input_fastas, &alphabet)?;

  // TODO: refactor to reduce duplication with `ancestral` as well as within the branches of this conditional
  if !dense {
    let graph: GraphAncestral = nwk_read_file(tree)?;

    #[allow(clippy::iter_on_single_items)]
    let partitions = [PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?, // TODO: allow other models
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }]
    .into_iter()
    .map(|p| Arc::new(RwLock::new(p)))
    .collect_vec();

    compress_sequences(&graph, &partitions, &aln)?;

    initial_guess_sparse(&graph, &partitions);
    let mut lh_prev = f64::MIN;
    for i in 0..*max_iter {
      let lh = run_marginal_sparse(&graph, &partitions)?;
      debug!("Iteration {}: likelihood {}", i + 1, float_to_significant_digits(lh, 7));
      if (lh - lh_prev).abs() < dp.abs() {
        break;
      }
      run_optimize_sparse(&graph, &partitions)?;
      lh_prev = lh;
    }

    write_graph(outdir, &graph)?;
  } else {
    let graph: GraphAncestral = nwk_read_file(tree)?;

    #[allow(clippy::iter_on_single_items)]
    let partitions = [PartitionMarginalDense {
      index: 0,
      gtr: jc69(JC69Params::default())?, // TODO: allow other models
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }]
    .into_iter()
    .map(|p| Arc::new(RwLock::new(p)))
    .collect_vec();

    let mut lh_prev = f64::MIN;
    for i in 0..*max_iter {
      let lh = run_marginal_dense(&graph, &partitions, &aln)?;

      // somehow, the initial guess makes it worse...
      // if i == 0 {
      //   initial_guess_dense(&graph, &partitions);
      // }

      debug!("Iteration {}: likelihood {}", i + 1, float_to_significant_digits(lh, 7));
      if (lh - lh_prev).abs() < dp.abs() {
        break;
      }
      run_optimize_dense(&graph, &partitions)?;
      lh_prev = lh;
    }

    write_graph(outdir, &graph)?;
  }

  Ok(())
}

fn write_graph<N, E, D>(outdir: impl AsRef<Path>, graph: &crate::graph::graph::Graph<N, E, D>) -> Result<(), Report>
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
