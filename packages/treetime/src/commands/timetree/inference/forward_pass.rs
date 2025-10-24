use crate::distribution::distribution_convolution::distribution_convolution;
use crate::distribution::distribution_division::distribution_division;
use crate::graph::breadth_first::GraphTraversalContinuation;
use crate::graph::graph::GraphNodeForward;
use crate::representation::graph_ancestral::{EdgeAncestral, GraphAncestral, NodeAncestral};
use eyre::Report;
use std::sync::Arc;

pub fn propagate_distributions_forward(graph: &GraphAncestral) -> Result<(), Report> {
  graph.par_iter_breadth_first_forward(|mut node| {
    propagate_distributions_forward_single_node(&mut node).unwrap();
    GraphTraversalContinuation::Continue
  });
  Ok(())
}

fn propagate_distributions_forward_single_node(
  node: &mut GraphNodeForward<NodeAncestral, EdgeAncestral, ()>,
) -> Result<(), Report> {
  set_likely_time(node);
  refine_distribution_from_parent(node)?;
  Ok(())
}

fn set_likely_time(node: &mut GraphNodeForward<NodeAncestral, EdgeAncestral, ()>) {
  let time = node
    .payload
    .time_distribution
    .as_ref()
    .and_then(|time_dist| time_dist.likely_time());

  node.payload.time = time;
}

fn refine_distribution_from_parent(
  node: &mut GraphNodeForward<NodeAncestral, EdgeAncestral, ()>,
) -> Result<(), Report> {
  if node.is_root {
    return Ok(());
  }

  if let Ok((parent, edge)) = node.get_exactly_one_parent() {
    let parent = parent.read_arc();
    let edge = edge.read_arc();

    if let (Some(parent_time_dist), Some(branch_dist), Some(subtree_dist)) =
      (&parent.time_distribution, &edge.branch_length_distribution, &node.payload.time_distribution)
    {
      let parent_except_subtree = if let Some(msg_to_parent) = &edge.msg_to_parent {
        distribution_division(parent_time_dist, msg_to_parent)?
      } else {
        parent_time_dist.as_ref().clone()
      };

      let dist_from_parent = distribution_convolution(&parent_except_subtree, branch_dist)?;
      let combined = distribution_convolution(&dist_from_parent, subtree_dist)?;
      node.payload.time_distribution = Some(Arc::new(combined));
    } else if let (Some(parent_time_dist), Some(branch_dist)) = (&parent.time_distribution, &edge.branch_length_distribution)
    {
      let dist_from_parent = distribution_convolution(parent_time_dist, branch_dist)?;
      node.payload.time_distribution = Some(Arc::new(dist_from_parent));
    }
  }

  Ok(())
}
