use crate::alphabet::alphabet::Alphabet;
use crate::commands::ancestral::fitch::compress_sequences;
use crate::commands::ancestral::marginal_dense::run_marginal_dense;
use crate::commands::ancestral::marginal_sparse::run_marginal_sparse;
use crate::commands::optimize::args::TreetimeOptimizeArgs;
use crate::commands::optimize::optimize_dense::run_optimize_dense;
use crate::commands::optimize::optimize_sparse::run_optimize_sparse;
use crate::graph::edge::GraphEdge;
use crate::graph::graph::Graph;
use crate::graph::node::GraphNode;
use crate::gtr::get_gtr::{get_gtr, get_gtr_dense};
use crate::io::fasta::read_many_fasta;
use crate::io::nex::{nex_write_file, NexWriteOptions};
use crate::io::nwk::{nwk_read_file, nwk_write_file, EdgeToNwk, NodeToNwk, NwkWriteOptions};
use crate::representation::graph_dense::DenseGraph;
use crate::representation::graph_sparse::SparseGraph;
use crate::representation::infer_dense::infer_dense;
use crate::representation::partitions_likelihood::{PartitionLikelihood, PartitionLikelihoodWithAln};
use crate::representation::partitions_parsimony::PartitionParsimonyWithAln;
use crate::utils::float_fmt::float_to_significant_digits;
use eyre::Report;
use itertools::Itertools;
use log::debug;
use serde::Serialize;
use std::path::Path;

use super::optimize_dense::initial_guess;

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
    let graph: SparseGraph = nwk_read_file(tree)?;
    let partitions = vec![PartitionParsimonyWithAln::new(alphabet.clone(), aln)?];
    let partitions = compress_sequences(&graph, partitions)?;

    let gtr = get_gtr(model_name, &alphabet, &graph)?;
    let partitions = partitions
          .into_iter()
          .map(|part| PartitionLikelihood::from_parsimony(gtr.clone(), part)) // FIXME: avoid cloning
          .collect_vec();

    let mut lh_prev = f64::MAX;
    for i in 0..*max_iter {
      let lh = run_marginal_sparse(&graph, &partitions)?;
      debug!("Iteration {}: likelihood {}", i + 1, float_to_significant_digits(lh, 7));
      if (lh_prev - lh).abs() < dp.abs() {
        break;
      }
      // run_optimize_sparse(&graph)?;
      lh_prev = lh;
    }

    write_graph(outdir, &graph)?;
  } else {
    let graph: DenseGraph = nwk_read_file(tree)?;
    let gtr = get_gtr_dense(model_name, &alphabet, &graph)?;

    let partitions_waln = vec![PartitionLikelihoodWithAln::new(gtr, alphabet, aln)?];
    let partitions = partitions_waln
      .iter()
      .map(|part| PartitionLikelihood::from(part.clone()))
      .collect_vec();
    let mut lh_prev = f64::MAX;
    for i in 0..*max_iter {
      //FIXME avoid assigning sequences to the graph in every iteration
      let lh = run_marginal_dense(&graph, partitions_waln.clone(), false)?; // FIXME: avoid cloning
                                                                            // if i == 0 {
                                                                            //   initial_guess(&graph, &partitions);
                                                                            // }
      debug!("Iteration {}: likelihood {}", i + 1, float_to_significant_digits(lh, 7));
      if (lh_prev - lh).abs() < dp.abs() {
        break;
      }
      run_optimize_dense(&graph, &partitions)?;
      lh_prev = lh;
    }

    write_graph(outdir, &graph)?;
  }

  Ok(())
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
