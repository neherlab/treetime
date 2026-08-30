use crate::ancestral::sample::SampleMode;
use crate::partition::traits::graph_log_lh;
use crate::partition::traits::{PartitionMarginalOps, PartitionMarginalPasses};
use eyre::Report;
use log::trace;
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_graph::edge::EdgeOptimizeOps;
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, Named};
use treetime_io::fasta::FastaRecord;
use treetime_primitives::{LogLh, Seq, seq};

pub fn initialize_marginal<N, E, P>(
  graph: &Graph<N, E, ()>,
  partitions: &[Arc<RwLock<P>>],
  aln: &[FastaRecord],
) -> Result<LogLh, Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
  P: PartitionMarginalOps<N, E> + ?Sized,
{
  for partition in partitions {
    partition.write_arc().attach_sequences(graph, aln)?;
  }
  update_marginal(graph, partitions)
}

pub fn update_marginal<N, E, P>(graph: &Graph<N, E, ()>, partitions: &[Arc<RwLock<P>>]) -> Result<LogLh, Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
  P: PartitionMarginalPasses<N, E> + ?Sized,
{
  marginal_backward(graph, partitions)?;
  let log_lh = graph_log_lh(graph, partitions)?;
  trace!("Marginal log likelihood (substitution): {}", log_lh.value());
  marginal_forward(graph, partitions)?;
  Ok(log_lh)
}

pub fn ancestral_reconstruction_marginal<N, E, P>(
  graph: &Graph<N, E, ()>,
  include_leaves: bool,
  impute: bool,
  partitions: &[Arc<RwLock<P>>],
  sample_mode: SampleMode,
  rng: &mut dyn rand::RngCore,
  mut visitor: impl FnMut(&N, &Seq) -> Result<(), Report>,
) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
  P: PartitionMarginalOps<N, E> + ?Sized,
{
  // Preorder traversal is sequential, so a single threaded RNG yields deterministic output under a
  // fixed seed: every node draws from the profile in a fixed traversal order.
  // Reconstruct every node so each partition's stored `seq.sequence` reflects the final, flag-aware
  // reconstruction that the node-data serializer reads back. `include_leaves` gates only whether tip
  // sequences are emitted to the visitor (the reconstructed FASTA), not whether they are computed:
  // `reconstruct_node_sequence` returns `None` for a suppressed tip after writing its `seq.sequence`.
  graph.iter_depth_first_preorder_forward(|node| {
    if partitions.is_empty() {
      if !include_leaves && node.is_leaf {
        return Ok(());
      }
      return visitor(&node.payload, &seq![]);
    }

    let reconstructed = {
      let mut partition = partitions[0].write_arc();
      partition.reconstruct_node_sequence(&node, include_leaves, impute, sample_mode, rng)
    };

    match reconstructed {
      Some(seq) => visitor(&node.payload, &seq),
      None => Ok(()),
    }
  })
}

pub fn marginal_backward<N, E, P>(graph: &Graph<N, E, ()>, partitions: &[Arc<RwLock<P>>]) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
  P: PartitionMarginalPasses<N, E> + ?Sized,
{
  for partition in partitions {
    let mut partition = partition.write_arc();
    partition.process_backward_pass(graph)?;
  }
  Ok(())
}

fn marginal_forward<N, E, P>(graph: &Graph<N, E, ()>, partitions: &[Arc<RwLock<P>>]) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
  P: PartitionMarginalPasses<N, E> + ?Sized,
{
  for partition in partitions {
    let mut partition = partition.write_arc();
    partition.process_forward_pass(graph)?;
  }
  Ok(())
}
