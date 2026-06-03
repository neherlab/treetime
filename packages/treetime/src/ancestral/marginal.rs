use crate::ancestral::sample::SampleMode;
use crate::partition::traits::graph_log_lh;
use crate::partition::traits::{PartitionMarginalOps, PartitionMarginalPasses};
use eyre::Report;
use log::trace;
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_graph::breadth_first::GraphTraversalContinuation;
use treetime_graph::edge::EdgeOptimizeOps;
use treetime_graph::graph::Graph;
use treetime_graph::graph_traverse::{GraphNodeBackward, GraphNodeForward};
use treetime_graph::node::{GraphNode, Named};
use treetime_io::fasta::FastaRecord;
use treetime_primitives::{Seq, seq};

pub fn initialize_marginal<N, E, P>(
  graph: &Graph<N, E, ()>,
  partitions: &[Arc<RwLock<P>>],
  aln: &[FastaRecord],
) -> Result<f64, Report>
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

pub fn update_marginal<N, E, P>(graph: &Graph<N, E, ()>, partitions: &[Arc<RwLock<P>>]) -> Result<f64, Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
  P: PartitionMarginalPasses<N, E> + ?Sized,
{
  marginal_backward(graph, partitions)?;
  let log_lh = graph_log_lh(graph, partitions)?;
  trace!("Marginal log likelihood (substitution): {log_lh}");
  marginal_forward(graph, partitions)?;
  Ok(log_lh)
}

pub fn ancestral_reconstruction_marginal<N, E, P>(
  graph: &Graph<N, E, ()>,
  include_leaves: bool,
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
  graph.iter_depth_first_preorder_forward(|node| {
    if !include_leaves && node.is_leaf {
      return Ok(());
    }

    let seq: Seq = if !partitions.is_empty() {
      let mut partition = partitions[0].write_arc();
      partition
        .reconstruct_node_sequence(&node, include_leaves, sample_mode, rng)
        .unwrap_or_else(|| {
          log::warn!("Missing reconstruction for node {:?}, using empty sequence", node.key);
          seq![]
        })
    } else {
      seq![]
    };

    visitor(&node.payload, &seq)
  })
}

pub fn marginal_backward<N, E, P>(graph: &Graph<N, E, ()>, partitions: &[Arc<RwLock<P>>]) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
  P: PartitionMarginalPasses<N, E> + ?Sized,
{
  graph.par_iter_breadth_first_backward(|node| {
    run_marginal_backward(partitions, &node)?;
    Ok(GraphTraversalContinuation::Continue)
  })
}

fn run_marginal_backward<N, E, P>(
  partitions: &[Arc<RwLock<P>>],
  node: &GraphNodeBackward<N, E, ()>,
) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
  P: PartitionMarginalPasses<N, E> + ?Sized,
{
  for partition in partitions {
    let mut partition = partition.write_arc();
    partition.process_node_backward(node)?;
  }
  Ok(())
}

fn marginal_forward<N, E, P>(graph: &Graph<N, E, ()>, partitions: &[Arc<RwLock<P>>]) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
  P: PartitionMarginalPasses<N, E> + ?Sized,
{
  graph.par_iter_breadth_first_forward(|node| {
    run_marginal_forward(graph, partitions, &node)?;
    Ok(GraphTraversalContinuation::Continue)
  })
}

fn run_marginal_forward<N, E, P>(
  graph: &Graph<N, E, ()>,
  partitions: &[Arc<RwLock<P>>],
  node: &GraphNodeForward<N, E, ()>,
) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
  P: PartitionMarginalPasses<N, E> + ?Sized,
{
  for partition in partitions {
    let mut partition = partition.write_arc();
    partition.process_node_forward(graph, node)?;
  }
  Ok(())
}
