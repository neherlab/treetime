use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::fitch::create_fitch_partition;
use crate::ancestral::gtr_inference::infer_gtr_fitch;
use crate::gtr::get_gtr::{GtrModelName, get_gtr_by_name, log_gtr};
use crate::gtr::gtr::GTR;
use crate::partition::algo::infer_dense::infer_dense;
use crate::partition::marginal_dense::PartitionMarginalDense;
use crate::partition::marginal_sparse::PartitionMarginalSparse;
use crate::seq::alignment::get_common_length;
use eyre::Report;
use treetime_graph::edge::{GraphEdge, HasBranchLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::NodeAncestralOps;
use treetime_io::fasta::FastaRecord;

pub enum MarginalPartition {
  Sparse(PartitionMarginalSparse),
  Dense(PartitionMarginalDense),
}

pub struct PartitionCreated {
  pub partition: MarginalPartition,
  pub gtr: GTR,
  pub model_name: GtrModelName,
}

/// Create a marginal partition from alignment data, consolidating the 3-way branch:
/// sparse, dense+infer GTR, dense+named GTR.
///
/// No file I/O. GTR JSON writing is the caller's responsibility.
pub fn create_marginal_partition<N, E>(
  graph: &Graph<N, E, ()>,
  index: usize,
  alphabet: Alphabet,
  sequences: &[FastaRecord],
  model_name: GtrModelName,
  dense: Option<bool>,
) -> Result<PartitionCreated, Report>
where
  N: NodeAncestralOps,
  E: GraphEdge + HasBranchLength,
{
  let dense = dense.unwrap_or_else(infer_dense);

  if !dense {
    let fitch = create_fitch_partition(graph, index, alphabet, sequences)?;
    let gtr = match model_name {
      GtrModelName::Infer => infer_gtr_fitch(&fitch, graph)?,
      _ => get_gtr_by_name(model_name)?,
    };
    log_gtr(&gtr, model_name);
    let partition = fitch.into_marginal_sparse(gtr.clone(), graph)?;
    Ok(PartitionCreated {
      partition: MarginalPartition::Sparse(partition),
      gtr,
      model_name,
    })
  } else if model_name == GtrModelName::Infer {
    let fitch = create_fitch_partition(graph, index, alphabet, sequences)?;
    let gtr = infer_gtr_fitch(&fitch, graph)?;
    log_gtr(&gtr, model_name);
    let partition = fitch.into_marginal_dense(gtr.clone());
    Ok(PartitionCreated {
      partition: MarginalPartition::Dense(partition),
      gtr,
      model_name,
    })
  } else {
    let length = get_common_length(sequences)?;
    let gtr = get_gtr_by_name(model_name)?;
    log_gtr(&gtr, model_name);
    let partition = PartitionMarginalDense::new(index, gtr.clone(), alphabet, length);
    Ok(PartitionCreated {
      partition: MarginalPartition::Dense(partition),
      gtr,
      model_name,
    })
  }
}
