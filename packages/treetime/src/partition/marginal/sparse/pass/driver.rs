use crate::partition::indexed_pass::IndexedPass;
use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
use crate::partition::marginal::sparse::pass::backward::process_node_backward_indexed;
use crate::partition::marginal::sparse::pass::forward::process_node_forward_indexed;
use eyre::Report;
use treetime_graph::edge::EdgeOptimizeOps;
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, Named};

pub fn process_backward_indexed<N, E>(
  partition: &mut PartitionMarginalSparse,
  graph: &Graph<N, E, ()>,
) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  let alphabet = partition.alphabet.clone();
  let gtr = partition.gtr.clone();
  let length = partition.length;
  let (nodes, edges) = (&mut partition.nodes, &mut partition.edges);
  let mut pass = IndexedPass::new(graph, nodes, edges, |key| {
    treetime_utils::make_internal_error!("Partition node {key} is missing before the sparse marginal pass")
  })?;
  let result = pass.try_for_each_backward(|dependencies, slot| {
    process_node_backward_indexed(graph, &alphabet, &gtr, length, dependencies, slot)
  });
  let (nodes, edges) = pass.into_maps()?;
  partition.nodes = nodes;
  partition.edges = edges;
  result
}

pub fn process_forward_indexed<N, E>(
  partition: &mut PartitionMarginalSparse,
  graph: &Graph<N, E, ()>,
) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  let alphabet = partition.alphabet.clone();
  let gtr = partition.gtr.clone();
  let length = partition.length;
  let root_sequence = partition.root_sequence.clone();
  let (nodes, edges) = (&mut partition.nodes, &mut partition.edges);
  let mut pass = IndexedPass::new(graph, nodes, edges, |key| {
    treetime_utils::make_internal_error!("Partition node {key} is missing before the sparse marginal pass")
  })?;
  let result = pass.try_for_each_forward(|dependencies, slot| {
    process_node_forward_indexed(graph, &alphabet, &gtr, length, &root_sequence, dependencies, slot)
  });
  let (nodes, edges) = pass.into_maps()?;
  partition.nodes = nodes;
  partition.edges = edges;
  result
}
