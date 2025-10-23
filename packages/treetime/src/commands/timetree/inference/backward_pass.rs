use crate::distribution::distribution::Distribution;
use crate::distribution::distribution_convolution::distribution_convolution;
use crate::graph::breadth_first::GraphTraversalContinuation;
use crate::graph::graph::GraphNodeBackward;
use crate::representation::graph_ancestral::{EdgeAncestral, GraphAncestral, NodeAncestral};
use eyre::Report;
use std::sync::Arc;

pub fn propagate_distributions_backward(graph: &GraphAncestral) -> Result<(), Report> {
  graph.par_iter_breadth_first_backward(|mut node| {
    propagate_distributions_backward_single_node(&mut node).unwrap();
    GraphTraversalContinuation::Continue
  });
  Ok(())
}

fn propagate_distributions_backward_single_node(
  node: &mut GraphNodeBackward<NodeAncestral, EdgeAncestral, ()>,
) -> Result<(), Report> {
  if node.is_leaf {
    return Ok(());
  }

  let mut result: Option<Distribution> = None;
  for (child, edge) in &node.children {
    let child = child.read_arc();
    let edge = edge.read_arc();

    if let (Some(branch_dist), Some(child_time_dist)) = (&edge.branch_length_distribution, &child.time_distribution) {
      let new = distribution_convolution(branch_dist, child_time_dist)?;
      result = Some(if let Some(result) = &result {
        distribution_convolution(result, &new)?
      } else {
        new
      });
    }
  }

  if let Some(dist) = result {
    node.payload.time_distribution = Some(Arc::new(dist));
  }

  Ok(())
}
