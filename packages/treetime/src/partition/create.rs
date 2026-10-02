use crate::alphabet::alphabet::Alphabet;
use crate::gtr::get_gtr::{GtrModelName, get_gtr_by_name, log_gtr};
use crate::gtr::gtr::GTR;
use crate::partition::algo::infer_dense::infer_dense;
use crate::partition::fitch::gtr_inference::infer_gtr_fitch;
use crate::partition::fitch::partition::PartitionFitch;
use crate::partition::fitch::passes::create_fitch_partition;
use crate::partition::marginal::dense::partition::PartitionMarginalDense;
use crate::partition::marginal::reconstruction::{DenseReconstruction, MarginalReconstruction, SparseReconstruction};
use crate::progress::LogSink;
use crate::seq::alignment::NodeSeqInput;
use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;

#[derive(Clone, Copy, Debug)]
pub(crate) enum Representation {
  Dense,
  Sparse,
}

impl Representation {
  pub(crate) fn resolve(dense: Option<bool>) -> Self {
    if dense.unwrap_or_else(infer_dense) {
      Self::Dense
    } else {
      Self::Sparse
    }
  }
}

pub(crate) fn build_marginal_partition(
  representation: Representation,
  model: GtrModelName,
  graph: &Graph,
  index: usize,
  alphabet: Alphabet,
  node_inputs: &BTreeMap<GraphNodeKey, NodeSeqInput>,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  log: &dyn LogSink,
) -> Result<MarginalReconstruction, Report> {
  Ok(match representation {
    Representation::Sparse => MarginalReconstruction::Sparse(build_sparse_reconstruction(
      model,
      graph,
      index,
      alphabet,
      node_inputs,
      branch_lengths,
      log,
    )?),
    Representation::Dense if model == GtrModelName::Infer => {
      let fitch = create_fitch_partition(graph, index, alphabet, node_inputs)?;
      let gtr = fitch_gtr(model, &fitch, graph, branch_lengths, log)?;
      let partition = fitch.into_marginal_dense(graph, node_inputs)?;
      MarginalReconstruction::Dense(DenseReconstruction::seeded(partition, gtr))
    },
    Representation::Dense => {
      let partition = PartitionMarginalDense::new(index, alphabet, graph, node_inputs)?;
      let gtr = named_gtr(model, log)?;
      MarginalReconstruction::Dense(DenseReconstruction::seeded(partition, gtr))
    },
  })
}

pub(crate) fn build_sparse_reconstruction(
  model: GtrModelName,
  graph: &Graph,
  index: usize,
  alphabet: Alphabet,
  node_inputs: &BTreeMap<GraphNodeKey, NodeSeqInput>,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  log: &dyn LogSink,
) -> Result<SparseReconstruction, Report> {
  let fitch = create_fitch_partition(graph, index, alphabet, node_inputs)?;
  let gtr = fitch_gtr(model, &fitch, graph, branch_lengths, log)?;
  let (partition, node_states) = fitch.into_marginal_sparse(graph)?;
  Ok(SparseReconstruction::seeded(partition, gtr, node_states))
}

fn fitch_gtr(
  model: GtrModelName,
  fitch: &PartitionFitch,
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  log: &dyn LogSink,
) -> Result<GTR, Report> {
  if model != GtrModelName::Infer {
    return named_gtr(model, log);
  }
  let gtr = infer_gtr_fitch(fitch, graph, branch_lengths, log)?;
  log_gtr(&gtr, model, log);
  Ok(gtr)
}

fn named_gtr(model: GtrModelName, log: &dyn LogSink) -> Result<GTR, Report> {
  let gtr = get_gtr_by_name(model)?;
  log_gtr(&gtr, model, log);
  Ok(gtr)
}
