use treetime_graph::breadth_first::GraphTraversalContinuation;
use treetime_graph::edge::EdgeOptimizeOps;
use treetime_graph::graph::Graph;
use treetime_graph::graph_traverse::{GraphNodeBackward, GraphNodeForward};
use treetime_graph::node::{GraphNode, Named};
use crate::io::fasta::FastaRecord;
use crate::representation::log_lh::HasLogLh;
use crate::representation::log_lh::graph_log_lh;
use crate::representation::partition_marginal::PartitionMarginalOps;
use crate::representation::seq::Seq;
use eyre::Report;
use log::debug;
use parking_lot::RwLock;
use std::sync::Arc;

/// Initialize partitions with sequence data and run marginal reconstruction.
///
/// Use this for the first marginal pass when sequences need to be attached.
/// For subsequent iterations, use `update_marginal` instead.
pub fn initialize_marginal<N, E, P>(
  graph: &Graph<N, E, ()>,
  partitions: &[Arc<RwLock<P>>],
  aln: &[FastaRecord],
) -> Result<f64, Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
  P: PartitionMarginalOps<N, E> + HasLogLh + ?Sized,
{
  for partition in partitions {
    partition.write_arc().attach_sequences(graph, aln)?;
  }
  update_marginal(graph, partitions)
}

/// Run marginal reconstruction without attaching sequences.
///
/// Use this for iterative refinement after sequences have been attached
/// via `initialize_marginal`.
pub fn update_marginal<N, E, P>(graph: &Graph<N, E, ()>, partitions: &[Arc<RwLock<P>>]) -> Result<f64, Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
  P: PartitionMarginalOps<N, E> + HasLogLh + ?Sized,
{
  marginal_backward(graph, partitions)?;
  let log_lh = graph_log_lh(graph, partitions)?;
  debug!("Log likelihood: {log_lh}");
  marginal_forward(graph, partitions)?;
  Ok(log_lh)
}

/// Ancestral sequence reconstruction
pub fn ancestral_reconstruction_marginal<N, E, P>(
  graph: &Graph<N, E, ()>,
  include_leaves: bool,
  partitions: &[Arc<RwLock<P>>],
  mut visitor: impl FnMut(&N, &Seq),
) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
  P: PartitionMarginalOps<N, E> + HasLogLh + ?Sized,
{
  graph.iter_depth_first_preorder_forward(|node| {
    if !include_leaves && node.is_leaf {
      return;
    }

    let seq: Seq = if !partitions.is_empty() {
      let mut partition = partitions[0].write_arc();
      partition
        .reconstruct_node_sequence(&node, include_leaves)
        .unwrap_or_else(|| crate::seq![])
    } else {
      crate::seq![]
    };

    visitor(&node.payload, &seq);
  });

  Ok(())
}

/// Backward pass: calculates ingroup profiles
fn marginal_backward<N, E, P>(graph: &Graph<N, E, ()>, partitions: &[Arc<RwLock<P>>]) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
  P: PartitionMarginalOps<N, E> + HasLogLh + ?Sized,
{
  graph.par_iter_breadth_first_backward(|node| {
    run_marginal_backward(partitions, &node).unwrap();
    GraphTraversalContinuation::Continue
  });
  Ok(())
}

fn run_marginal_backward<N, E, P>(
  partitions: &[Arc<RwLock<P>>],
  node: &GraphNodeBackward<N, E, ()>,
) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
  P: PartitionMarginalOps<N, E> + HasLogLh + ?Sized,
{
  for partition in partitions {
    let mut partition = partition.write_arc();
    partition.process_node_backward(node)?;
  }
  Ok(())
}

/// Forward pass: calculates outgroup profiles
fn marginal_forward<N, E, P>(graph: &Graph<N, E, ()>, partitions: &[Arc<RwLock<P>>]) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
  P: PartitionMarginalOps<N, E> + HasLogLh + ?Sized,
{
  graph.par_iter_breadth_first_forward(|node| {
    run_marginal_forward(graph, partitions, &node).unwrap();
    GraphTraversalContinuation::Continue
  });
  Ok(())
}

fn run_marginal_forward<N, E, P>(
  graph: &Graph<N, E, ()>,
  partitions: &[Arc<RwLock<P>>],
  node: &GraphNodeForward<N, E, ()>,
) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
  P: PartitionMarginalOps<N, E> + HasLogLh + ?Sized,
{
  for partition in partitions {
    let mut partition = partition.write_arc();
    partition.process_node_forward(graph, node)?;
  }
  Ok(())
}
