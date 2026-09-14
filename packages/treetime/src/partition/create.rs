use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::fitch::create_fitch_partition;
use crate::ancestral::gtr_inference::infer_gtr_fitch;
use crate::gtr::get_gtr::{GtrModelName, get_gtr_by_name, log_gtr};
use crate::partition::algo::infer_dense::infer_dense;
use crate::partition::marginal::dense::partition::PartitionMarginalDense;
use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
use crate::partition::storage::sparse::SparseNodeState;
use crate::seq::alignment::get_common_length_of_node_inputs;
use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::nwk::NwkFastaNodeInput;

pub enum MarginalPartition {
  /// A sparse partition together with the seed node-state map derived from the Fitch handoff. The
  /// caller owns and threads the node states; the partition holds only the durable inputs.
  Sparse(PartitionMarginalSparse, BTreeMap<GraphNodeKey, SparseNodeState>),
  Dense(PartitionMarginalDense),
}

pub struct PartitionCreated {
  pub partition: MarginalPartition,
  pub model_name: GtrModelName,
}

/// Create a marginal partition from alignment data, consolidating the 3-way branch:
/// sparse, dense+infer GTR, dense+named GTR.
///
/// No file I/O. GTR JSON writing is the caller's responsibility.
pub fn create_marginal_partition(
  graph: &Graph,
  index: usize,
  alphabet: Alphabet,
  node_inputs: &BTreeMap<GraphNodeKey, NwkFastaNodeInput>,
  model_name: GtrModelName,
  dense: Option<bool>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<PartitionCreated, Report> {
  let dense = dense.unwrap_or_else(infer_dense);

  if !dense {
    let fitch = create_fitch_partition(graph, index, alphabet, node_inputs)?;
    let gtr = match model_name {
      GtrModelName::Infer => infer_gtr_fitch(&fitch, graph, branch_lengths)?,
      _ => get_gtr_by_name(model_name)?,
    };
    log_gtr(&gtr, model_name);
    let (partition, node_states) = fitch.into_marginal_sparse(gtr, graph)?;
    Ok(PartitionCreated {
      partition: MarginalPartition::Sparse(partition, node_states),
      model_name,
    })
  } else if model_name == GtrModelName::Infer {
    let fitch = create_fitch_partition(graph, index, alphabet, node_inputs)?;
    let gtr = infer_gtr_fitch(&fitch, graph, branch_lengths)?;
    log_gtr(&gtr, model_name);
    let partition = fitch.into_marginal_dense(gtr);
    Ok(PartitionCreated {
      partition: MarginalPartition::Dense(partition),
      model_name,
    })
  } else {
    let length = get_common_length_of_node_inputs(node_inputs)?;
    let gtr = get_gtr_by_name(model_name)?;
    log_gtr(&gtr, model_name);
    let partition = PartitionMarginalDense::new(index, gtr, alphabet, length);
    Ok(PartitionCreated {
      partition: MarginalPartition::Dense(partition),
      model_name,
    })
  }
}
