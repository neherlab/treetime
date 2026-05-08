use crate::representation::partition::traits::HasLogLh;
use crate::representation::partition::traits::PartitionMarginalOps;
use crate::representation::partition::traits::graph_log_lh;
use eyre::Report;
use log::trace;
use parking_lot::{Mutex, RwLock};
use std::sync::Arc;
use treetime_graph::breadth_first::GraphTraversalContinuation;
use treetime_graph::edge::EdgeOptimizeOps;
use treetime_graph::graph::Graph;
use treetime_graph::graph_traverse::{GraphNodeBackward, GraphNodeForward};
use treetime_graph::node::{GraphNode, Named};
use treetime_io::fasta::FastaRecord;
use treetime_primitives::{Seq, seq};
use treetime_utils::sync::mutex::extract_parallel_error;

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
  trace!("Marginal log likelihood (substitution): {log_lh}");
  marginal_forward(graph, partitions)?;
  Ok(log_lh)
}

/// Ancestral sequence reconstruction
pub fn ancestral_reconstruction_marginal<N, E, P>(
  graph: &Graph<N, E, ()>,
  include_leaves: bool,
  partitions: &[Arc<RwLock<P>>],
  mut visitor: impl FnMut(&N, &Seq) -> Result<(), Report>,
) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
  P: PartitionMarginalOps<N, E> + HasLogLh + ?Sized,
{
  graph.try_iter_depth_first_preorder_forward(|node| {
    if !include_leaves && node.is_leaf {
      return Ok(());
    }

    let seq: Seq = if !partitions.is_empty() {
      let mut partition = partitions[0].write_arc();
      partition
        .reconstruct_node_sequence(&node, include_leaves)
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

/// Backward pass: calculates ingroup profiles
fn marginal_backward<N, E, P>(graph: &Graph<N, E, ()>, partitions: &[Arc<RwLock<P>>]) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
  P: PartitionMarginalOps<N, E> + HasLogLh + ?Sized,
{
  let error: Arc<Mutex<Option<Report>>> = Arc::new(Mutex::new(None));
  graph.par_iter_breadth_first_backward(|node| {
    if let Err(e) = run_marginal_backward(partitions, &node) {
      let mut guard = error.lock();
      if guard.is_none() {
        *guard = Some(e);
      }
      return GraphTraversalContinuation::Stop;
    }
    GraphTraversalContinuation::Continue
  });
  extract_parallel_error(error)
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
  let error: Arc<Mutex<Option<Report>>> = Arc::new(Mutex::new(None));
  graph.par_iter_breadth_first_forward(|node| {
    if let Err(e) = run_marginal_forward(graph, partitions, &node) {
      let mut guard = error.lock();
      if guard.is_none() {
        *guard = Some(e);
      }
      return GraphTraversalContinuation::Stop;
    }
    GraphTraversalContinuation::Continue
  });
  extract_parallel_error(error)
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
