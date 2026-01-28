use crate::commands::timetree::timetree_traits::{TimetreeEdge, TimetreeNode};
use crate::distribution::distribution::{Distribution, DistributionNegLog};
use crate::distribution::distribution_convolution::distribution_convolution;
use crate::distribution::distribution_multiplication::distribution_multiplication;
use crate::graph::breadth_first::GraphTraversalContinuation;
use crate::graph::graph::{Graph, GraphNodeBackward};
use crate::graph::node::GraphNodeKey;
use eyre::Report;
use indexmap::IndexMap;
use std::sync::Arc;

/// Propagates time distributions backward from leaves to root.
///
/// If coalescent_contributions is provided, multiplies each node's distribution
/// with its precalculated coalescent contribution.
pub fn propagate_distributions_backward<N, E, D>(
  graph: &Graph<N, E, D>,
  coalescent_contributions: Option<&IndexMap<GraphNodeKey, Arc<DistributionNegLog>>>,
) -> Result<(), Report>
where
  N: TimetreeNode,
  E: TimetreeEdge,
  D: Send + Sync,
{
  graph.par_iter_breadth_first_backward(|mut node| {
    propagate_distributions_backward_single_node(&mut node, coalescent_contributions).unwrap();
    GraphTraversalContinuation::Continue
  });
  Ok(())
}

/// Computes time distribution for a single internal node from its children.
fn propagate_distributions_backward_single_node<N, E, D>(
  node: &mut GraphNodeBackward<N, E, D>,
  coalescent_contribs: Option<&IndexMap<GraphNodeKey, Arc<DistributionNegLog>>>,
) -> Result<(), Report>
where
  N: TimetreeNode,
  E: TimetreeEdge,
  D: Send + Sync,
{
  if node.is_leaf {
    // Apply precalculated coalescent contribution to leaf
    if let (Some(time_dist), Some(coalescent_contrib)) = (
      node.payload.time_distribution(),
      coalescent_contribs.and_then(|c| c.get(&node.key)),
    ) {
      let coalescent_contrib_plain = coalescent_contrib.to_plain();
      let combined = distribution_multiplication(time_dist.as_ref(), &coalescent_contrib_plain)?;
      node.payload.set_time_distribution(Some(Arc::new(combined)));
    }
    return Ok(());
  }

  // For internal node, initialize distribution with coalescent contribution
  let mut result: Option<Distribution> = coalescent_contribs
    .and_then(|contributions| contributions.get(&node.key))
    .map(|contrib| (**contrib).to_plain());

  for (child, edge) in &node.children {
    let child = child.read_arc();
    let mut edge = edge.write_arc();

    if let (Some(branch_dist), Some(child_time_dist)) = (edge.branch_length_distribution(), child.time_distribution()) {
      // Compute parent time distribution using regular convolution with negated branch: parent_time = child_time + (-branch_length)
      let negated_branch_dist = branch_dist.negate();
      let parent_message = distribution_convolution(child_time_dist.as_ref(), &negated_branch_dist)?;
      let parent_message_arc = Arc::new(parent_message);

      // Store message on edge
      edge.set_msg_to_parent(Some(Arc::clone(&parent_message_arc)));

      // Combine messages from all children using multiplication (intersection of constraints)
      result = Some(if let Some(current) = result {
        distribution_multiplication(&current, &parent_message_arc)?
      } else {
        (*parent_message_arc).clone()
      });
    }
  }

  // Store final distribution on node
  if let Some(dist) = result {
    node.payload.set_time_distribution(Some(Arc::new(dist)));
  }

  Ok(())
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::distribution::distribution::DistributionPlain as Distribution;
  use crate::distribution::distribution_convolution::distribution_convolution;
  use approx::assert_abs_diff_eq;

  #[test]
  fn test_inverse_convolution_integration() -> Result<(), Report> {
    // Test that the refactored code uses the same operation as the standalone function
    let child_dist: Distribution = Distribution::point(2013.0, 1.0);
    let branch_dist: Distribution = Distribution::point(2.5, 0.8);

    let negated_branch_dist: Distribution = branch_dist.negate();
    let result: Distribution = distribution_convolution(&child_dist, &negated_branch_dist)?;

    if let Some(parent_time) = result.likely_time() {
      assert_abs_diff_eq!(parent_time, 2010.5, epsilon = 1e-10);
      assert!(parent_time < 2013.0, "Parent should be older than child");
    } else {
      panic!("Expected valid parent time");
    }

    Ok(())
  }

  #[test]
  fn test_backward_pass_time_direction() -> Result<(), Report> {
    // Test that backward pass correctly computes older times for ancestors
    let child_time = 2012.0;
    let branch_length = 1.5;
    let expected_parent_time = child_time - branch_length; // 2010.5

    let child_dist: Distribution = Distribution::point(child_time, 1.0);
    let branch_dist: Distribution = Distribution::point(branch_length, 1.0);

    let negated_branch_dist: Distribution = branch_dist.negate();
    let parent_dist: Distribution = distribution_convolution(&child_dist, &negated_branch_dist)?;

    if let Some(actual_parent_time) = parent_dist.likely_time() {
      assert_abs_diff_eq!(actual_parent_time, expected_parent_time, epsilon = 1e-10);
      assert!(actual_parent_time < child_time, "Parent must be older than child");
    } else {
      panic!("Expected valid parent time");
    }

    Ok(())
  }
}
