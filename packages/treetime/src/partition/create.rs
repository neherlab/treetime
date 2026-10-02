use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::fitch::create_fitch_partition;
use crate::ancestral::gtr_inference::infer_gtr_fitch;
use crate::ancestral::plan::Representation;
use crate::gtr::get_gtr::{GtrModelName, get_gtr_by_name, log_gtr};
use crate::gtr::gtr::GTR;
use crate::partition::fitch::partition::PartitionFitch;
use crate::partition::marginal::dense::partition::PartitionMarginalDense;
use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
use crate::partition::storage::sparse::SparseNodeState;
use crate::progress::LogSink;
use crate::seq::alignment::NodeSeqInput;
use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;

pub(crate) fn build_marginal_partition(
  representation: Representation,
  model: GtrModelName,
  graph: &Graph,
  index: usize,
  alphabet: Alphabet,
  node_inputs: &BTreeMap<GraphNodeKey, NodeSeqInput>,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  log: &dyn LogSink,
) -> Result<(MarginalPartition, GTR), Report> {
  match representation {
    Representation::Sparse => {
      let (partition, node_states, gtr) =
        build_sparse_partition(model, graph, index, alphabet, node_inputs, branch_lengths, log)?;
      Ok((MarginalPartition::Sparse(partition, node_states), gtr))
    },
    Representation::Dense => {
      let (partition, gtr) = build_dense_partition(model, graph, index, alphabet, node_inputs, branch_lengths, log)?;
      Ok((MarginalPartition::Dense(partition), gtr))
    },
  }
}

pub(crate) fn build_sparse_partition(
  model: GtrModelName,
  graph: &Graph,
  index: usize,
  alphabet: Alphabet,
  node_inputs: &BTreeMap<GraphNodeKey, NodeSeqInput>,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  log: &dyn LogSink,
) -> Result<(PartitionMarginalSparse, BTreeMap<GraphNodeKey, SparseNodeState>, GTR), Report> {
  let fitch = create_fitch_partition(graph, index, alphabet, node_inputs)?;
  let gtr = build_gtr(model, Some(&fitch), graph, branch_lengths, log)?;
  let (partition, node_states) = fitch.into_marginal_sparse(graph)?;
  Ok((partition, node_states, gtr))
}

fn build_dense_partition(
  model: GtrModelName,
  graph: &Graph,
  index: usize,
  alphabet: Alphabet,
  node_inputs: &BTreeMap<GraphNodeKey, NodeSeqInput>,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  log: &dyn LogSink,
) -> Result<(PartitionMarginalDense, GTR), Report> {
  if model == GtrModelName::Infer {
    let fitch = create_fitch_partition(graph, index, alphabet, node_inputs)?;
    let gtr = build_gtr(model, Some(&fitch), graph, branch_lengths, log)?;
    Ok((fitch.into_marginal_dense(graph, node_inputs)?, gtr))
  } else {
    let partition = PartitionMarginalDense::new(index, alphabet, graph, node_inputs)?;
    let gtr = build_gtr(model, None, graph, branch_lengths, log)?;
    Ok((partition, gtr))
  }
}

fn build_gtr(
  model: GtrModelName,
  fitch: Option<&PartitionFitch>,
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  log: &dyn LogSink,
) -> Result<GTR, Report> {
  let gtr = match fitch {
    Some(fitch) if model == GtrModelName::Infer => infer_gtr_fitch(fitch, graph, branch_lengths, log)?,
    _ => get_gtr_by_name(model)?,
  };
  log_gtr(&gtr, model, log);
  Ok(gtr)
}

pub(crate) enum MarginalPartition {
  Sparse(PartitionMarginalSparse, BTreeMap<GraphNodeKey, SparseNodeState>),
  Dense(PartitionMarginalDense),
}
