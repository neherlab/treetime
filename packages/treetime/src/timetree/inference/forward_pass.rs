use crate::payload::traits::{TimetreeEdge, TimetreeNode};
use eyre::Report;
use std::sync::Arc;
use treetime_distribution::distribution_convolution;
use treetime_distribution::distribution_division;
use treetime_distribution::distribution_multiplication;
use treetime_graph::breadth_first::GraphTraversalContinuation;
use treetime_graph::graph::Graph;
use treetime_graph::graph_traverse::GraphNodeForward;

pub fn propagate_distributions_forward<N, E, D>(graph: &Graph<N, E, D>) -> Result<(), Report>
where
  N: TimetreeNode,
  E: TimetreeEdge,
  D: Send + Sync,
{
  graph.par_iter_breadth_first_forward(|mut node| {
    propagate_distributions_forward_single_node(&mut node)?;
    Ok(GraphTraversalContinuation::Continue)
  })
}

fn propagate_distributions_forward_single_node<N, E, D>(node: &mut GraphNodeForward<N, E, D>) -> Result<(), Report>
where
  N: TimetreeNode,
  E: TimetreeEdge,
  D: Send + Sync,
{
  refine_distribution_from_parent(node)?;
  set_likely_time(node);
  Ok(())
}

fn set_likely_time<N, E, D>(node: &mut GraphNodeForward<N, E, D>)
where
  N: TimetreeNode,
  E: TimetreeEdge,
  D: Send + Sync,
{
  let time = node
    .payload
    .time_distribution()
    .as_ref()
    .and_then(|time_dist| time_dist.likely_time());

  if let Some(time) = time {
    node.payload.set_time(Some(time));
  }
}

fn refine_distribution_from_parent<N, E, D>(node: &mut GraphNodeForward<N, E, D>) -> Result<(), Report>
where
  N: TimetreeNode,
  E: TimetreeEdge,
  D: Send + Sync,
{
  if node.is_root {
    return Ok(());
  }

  // Do not overwrite leaf time_distribution (date constraint)
  if node.is_leaf {
    return Ok(());
  }

  if let Ok((parent, edge)) = node.get_exactly_one_parent() {
    let parent = parent.read_arc();
    let edge = edge.read_arc();

    if let (Some(parent_time_dist), Some(branch_dist), Some(subtree_dist)) = (
      parent.time_distribution(),
      edge.branch_length_distribution(),
      node.payload.time_distribution(),
    ) {
      let parent_except_subtree = if let Some(msg_to_parent) = edge.msg_to_parent() {
        distribution_division(parent_time_dist, msg_to_parent)?
      } else {
        parent_time_dist.as_ref().clone()
      };

      let dist_from_parent = distribution_convolution(&parent_except_subtree, branch_dist)?;
      // Normalize to prevent numerical underflow: the backward pass stores normalized
      // distributions (max=1.0), and the convolution/division can produce arbitrary scales.
      // Without normalization, values accumulate downward across tree depth.
      let combined = distribution_multiplication(&dist_from_parent, subtree_dist)?.normalize();
      node.payload.set_time_distribution(Some(Arc::new(combined)));
    } else if let (Some(parent_time_dist), Some(branch_dist)) =
      (parent.time_distribution(), edge.branch_length_distribution())
    {
      let dist_from_parent = distribution_convolution(parent_time_dist, branch_dist)?.normalize();
      node.payload.set_time_distribution(Some(Arc::new(dist_from_parent)));
    }
  }

  Ok(())
}
