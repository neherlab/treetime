use crate::commands::timetree::inference::utils::convolve_safe;
use crate::graph::breadth_first::GraphTraversalContinuation;
use crate::graph::graph::GraphNodeForward;
use crate::representation::graph_ancestral::{EdgeAncestral, GraphAncestral, NodeAncestral};
use eyre::Report;
use std::sync::Arc;

pub fn propagate_distributions_forward(graph: &GraphAncestral) -> Result<(), Report> {
  graph.par_iter_breadth_first_forward(|mut node| {
    set_likely_time(&mut node);
    refine_distribution_from_parent(&mut node);
    GraphTraversalContinuation::Continue
  });
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

fn refine_distribution_from_parent(node: &mut GraphNodeForward<NodeAncestral, EdgeAncestral, ()>) {
  if node.is_root {
    return;
  }

  if let Ok((parent, edge)) = node.get_exactly_one_parent() {
    let parent = parent.read_arc();
    let edge = edge.read_arc();

    if let (Some(parent_time_dist), Some(branch_dist)) = (&parent.time_distribution, &edge.branch_length_distribution) {
      let dist_from_parent = convolve_safe(parent_time_dist, branch_dist);

      let combined = if let Some(subtree_dist) = &node.payload.time_distribution {
        convolve_safe(&dist_from_parent, subtree_dist)
      } else {
        dist_from_parent
      };

      node.payload.time_distribution = Some(Arc::new(combined));
    }
  }
}
