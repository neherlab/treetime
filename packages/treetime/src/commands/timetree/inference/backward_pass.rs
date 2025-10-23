use crate::commands::timetree::inference::utils::{convolve_safe, multiply_distributions};
use crate::distribution::distribution::Distribution;
use crate::graph::breadth_first::GraphTraversalContinuation;
use crate::representation::graph_ancestral::GraphAncestral;
use eyre::Report;
use std::sync::Arc;

pub fn propagate_distributions_backward(graph: &GraphAncestral) -> Result<(), Report> {
  graph.par_iter_breadth_first_backward(|mut node| {
    if !node.is_leaf {
      let mut result: Option<Distribution> = None;

      for (child, edge) in &node.children {
        let child = child.read_arc();
        let edge = edge.read_arc();

        if let (Some(branch_dist), Some(child_time_dist)) = (&edge.branch_length_distribution, &child.time_distribution)
        {
          let convolved = convolve_safe(branch_dist, child_time_dist);
          result = Some(multiply_distributions(result, convolved));
        }
      }

      if let Some(dist) = result {
        node.payload.time_distribution = Some(Arc::new(dist));
      }
    }
    GraphTraversalContinuation::Continue
  });
  Ok(())
}
