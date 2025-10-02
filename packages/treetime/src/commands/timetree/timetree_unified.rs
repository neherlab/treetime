use crate::graph::breadth_first::GraphTraversalContinuation;
use crate::graph::graph::{GraphNodeBackward, GraphNodeForward};
use crate::representation::graph_ancestral::{EdgeAncestral, GraphAncestral, NodeAncestral};
use crate::representation::partition_timetree::PartitionTimetreeOps;
use eyre::Report;
use parking_lot::RwLock;
use std::sync::Arc;

pub fn run_timetree(
  graph: &GraphAncestral,
  partitions: &[Arc<RwLock<dyn PartitionTimetreeOps>>],
) -> Result<(), Report> {
  timetree_backward(graph, partitions)?;
  timetree_forward(graph, partitions)?;
  Ok(())
}

fn timetree_backward(
  graph: &GraphAncestral,
  partitions: &[Arc<RwLock<dyn PartitionTimetreeOps>>],
) -> Result<(), Report> {
  graph.par_iter_breadth_first_backward(|node| {
    run_timetree_backward(partitions, &node).unwrap();
    GraphTraversalContinuation::Continue
  });
  Ok(())
}

fn run_timetree_backward(
  partitions: &[Arc<RwLock<dyn PartitionTimetreeOps>>],
  node: &GraphNodeBackward<NodeAncestral, EdgeAncestral, ()>,
) -> Result<(), Report> {
  for partition in partitions {
    let mut partition = partition.write_arc();
    partition.process_node_backward(node)?;
  }
  Ok(())
}

fn timetree_forward(
  graph: &GraphAncestral,
  partitions: &[Arc<RwLock<dyn PartitionTimetreeOps>>],
) -> Result<(), Report> {
  graph.par_iter_breadth_first_forward(|node| {
    run_timetree_forward(graph, partitions, &node).unwrap();
    GraphTraversalContinuation::Continue
  });
  Ok(())
}

fn run_timetree_forward(
  graph: &GraphAncestral,
  partitions: &[Arc<RwLock<dyn PartitionTimetreeOps>>],
  node: &GraphNodeForward<NodeAncestral, EdgeAncestral, ()>,
) -> Result<(), Report> {
  for partition in partitions {
    let mut partition = partition.write_arc();
    partition.process_node_forward(graph, node)?;
  }
  Ok(())
}
