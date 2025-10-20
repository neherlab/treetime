use crate::graph::breadth_first::GraphTraversalContinuation;
use crate::graph::graph::{GraphNodeBackward, GraphNodeForward};
use crate::representation::edge_timetree::EdgeTimetree;
use crate::representation::node_timetree::NodeTimetree;
use crate::representation::partition_timetree::{GraphTimetree, PartitionTimetreeOps};
use eyre::Report;
use parking_lot::RwLock;
use std::sync::Arc;

pub fn run_timetree<P>(
  graph: &GraphTimetree,
  partitions: &[Arc<RwLock<P>>],
) -> Result<(), Report>
where
  P: PartitionTimetreeOps<NodeTimetree, EdgeTimetree> + ?Sized,
{
  timetree_backward(graph, partitions)?;
  timetree_forward(graph, partitions)?;
  Ok(())
}

fn timetree_backward<P>(
  graph: &GraphTimetree,
  partitions: &[Arc<RwLock<P>>],
) -> Result<(), Report>
where
  P: PartitionTimetreeOps<NodeTimetree, EdgeTimetree> + ?Sized,
{
  graph.par_iter_breadth_first_backward(|node| {
    run_timetree_backward(partitions, &node).unwrap();
    GraphTraversalContinuation::Continue
  });
  Ok(())
}

fn run_timetree_backward<P>(
  partitions: &[Arc<RwLock<P>>],
  node: &GraphNodeBackward<NodeTimetree, EdgeTimetree, ()>,
) -> Result<(), Report>
where
  P: PartitionTimetreeOps<NodeTimetree, EdgeTimetree> + ?Sized,
{
  for partition in partitions {
    let mut partition = partition.write_arc();
    partition.process_node_backward(node)?;
  }
  Ok(())
}

fn timetree_forward<P>(
  graph: &GraphTimetree,
  partitions: &[Arc<RwLock<P>>],
) -> Result<(), Report>
where
  P: PartitionTimetreeOps<NodeTimetree, EdgeTimetree> + ?Sized,
{
  graph.par_iter_breadth_first_forward(|node| {
    run_timetree_forward(graph, partitions, &node).unwrap();
    GraphTraversalContinuation::Continue
  });
  Ok(())
}

fn run_timetree_forward<P>(
  graph: &GraphTimetree,
  partitions: &[Arc<RwLock<P>>],
  node: &GraphNodeForward<NodeTimetree, EdgeTimetree, ()>,
) -> Result<(), Report>
where
  P: PartitionTimetreeOps<NodeTimetree, EdgeTimetree> + ?Sized,
{
  for partition in partitions {
    let mut partition = partition.write_arc();
    partition.process_node_forward(graph, node)?;
  }
  Ok(())
}
