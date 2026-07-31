use crate::partition::indexed_pass::{IndexedPassDependencies, IndexedPassSlot, with_indexed_graph_payloads};
use crate::payload::traits::{TimetreeEdge, TimetreeNode};
use eyre::Report;
use log::warn;
use std::sync::Arc;
use treetime_distribution::BoundaryBehavior;
use treetime_distribution::distribution_convolution;
use treetime_distribution::distribution_division;
use treetime_distribution::distribution_multiplication;
use treetime_graph::edge::GraphEdge;
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, Named};

pub fn propagate_distributions_forward<N, E, D>(graph: &Graph<N, E, D>) -> Result<(), Report>
where
  N: GraphNode + Named + TimetreeNode + Default,
  E: GraphEdge + TimetreeEdge + Default,
  D: Send + Sync,
{
  with_indexed_graph_payloads(graph, |pass| {
    pass.try_for_each_forward(|dependencies, slot| propagate_distributions_forward_slot(graph, dependencies, slot))
  })
}

fn propagate_distributions_forward_slot<N, E, D>(
  graph: &Graph<N, E, D>,
  dependencies: &IndexedPassDependencies<N, E>,
  slot: &mut IndexedPassSlot<N, E>,
) -> Result<(), Report>
where
  N: GraphNode + Named + TimetreeNode,
  E: GraphEdge + TimetreeEdge,
  D: Send + Sync,
{
  refine_distribution_from_parent(graph, dependencies, slot)?;

  let is_leaf = graph.is_leaf(slot.key);

  // Independent marginal modes can invert even though the forward pass commits each
  // parent before visiting its children. Current behavior projects non-leaf internal
  // point estimates to the committed parent time; leaves retain their observed dates.
  // This changes the point estimate without recomputing the posterior. The statistical
  // contract remains open in kb/issues/M-timetree-marginal-node-times-can-violate-topology.md.
  let parent_time = (!is_leaf).then(|| parent_time(dependencies, slot)).flatten();

  // An internal node whose time distribution is empty (or degenerate) yields no likely
  // time, so no date is assigned. This signals irreconcilable inference constraints at the
  // node (e.g. disjoint child messages). Surface it instead of dropping the date silently.
  if set_likely_time(&mut slot.node, parent_time).is_none() && !is_leaf {
    let name = slot.node.name();
    let name = name.as_ref().map_or("<unnamed>", |name| name.as_ref());
    warn!(
      "Timetree forward pass: internal node '{name}' has an empty time distribution; no date was \
       assigned. The inference constraints at this node are likely irreconcilable (e.g. conflicting \
       or disjoint date constraints from its subtree)."
    );
  }
  Ok(())
}

/// Committed time of the slot's parent, if it has one and it is set.
fn parent_time<N, E>(dependencies: &IndexedPassDependencies<N, E>, slot: &IndexedPassSlot<N, E>) -> Option<f64>
where
  N: GraphNode + TimetreeNode,
  E: GraphEdge + TimetreeEdge,
{
  let parent_key = slot.parent_key?;
  dependencies.node(parent_key).time()
}

/// Assign the node's committed time from the peak of its time distribution, clamped to be no
/// earlier than the parent's committed time. Returns the assigned time, or `None` when the time
/// distribution is empty or degenerate and no time could be determined (the node is left undated).
fn set_likely_time(node: &mut impl TimetreeNode, parent_time: Option<f64>) -> Option<f64> {
  let time = node
    .time_distribution()
    .as_ref()
    .and_then(|time_dist| time_dist.likely_time())?;

  let time = parent_time.map_or(time, |parent_time| time.max(parent_time));
  node.set_time(Some(time));
  Some(time)
}

fn refine_distribution_from_parent<N, E, D>(
  graph: &Graph<N, E, D>,
  dependencies: &IndexedPassDependencies<N, E>,
  slot: &mut IndexedPassSlot<N, E>,
) -> Result<(), Report>
where
  N: GraphNode + TimetreeNode,
  E: GraphEdge + TimetreeEdge,
  D: Send + Sync,
{
  if slot.parent_key.is_none() {
    return Ok(());
  }

  // Do not overwrite leaf time_distribution (date constraint)
  if graph.is_leaf(slot.key) {
    return Ok(());
  }

  if let Some(parent_key) = slot.parent_key {
    let parent = dependencies.node(parent_key);
    let (_, edge) = slot
      .parent_edge
      .as_ref()
      .expect("Non-root indexed node must own its parent edge");

    if let (Some(parent_time_dist), Some(branch_dist), Some(subtree_dist)) = (
      parent.time_distribution(),
      edge.branch_length_distribution(),
      slot.node.time_distribution(),
    ) {
      let parent_except_subtree = if let Some(msg_to_parent) = edge.msg_to_parent() {
        distribution_division(parent_time_dist, msg_to_parent)?
      } else {
        parent_time_dist.as_ref().clone()
      };

      // Tail policy for the forward message (kb/decisions/distribution-tails-and-arithmetic.md).
      // The parent's time is a hard lower bound (left tail Zero); there is no upper bound from
      // the parent side on how far in the future the node could be (right tail Constant).
      let dist_from_parent = distribution_convolution(&parent_except_subtree, branch_dist)?
        .with_left_extrap(BoundaryBehavior::Zero)?
        .with_right_extrap(BoundaryBehavior::Constant)?;
      // Normalize to prevent numerical underflow: the backward pass stores normalized
      // distributions (max=1.0), and the convolution/division can produce arbitrary scales.
      // Without normalization, values accumulate downward across tree depth.
      let combined = distribution_multiplication(&dist_from_parent, subtree_dist)?.normalize();
      slot.node.set_time_distribution(Some(Arc::new(combined)));
    } else if let (Some(parent_time_dist), Some(branch_dist)) =
      (parent.time_distribution(), edge.branch_length_distribution())
    {
      // Same forward tail policy as the two-parent-message branch above. normalize() rebuilds
      // the grid, so the tail policy is applied afterwards, not before.
      let dist_from_parent = distribution_convolution(parent_time_dist, branch_dist)?
        .normalize()
        .with_left_extrap(BoundaryBehavior::Zero)?
        .with_right_extrap(BoundaryBehavior::Constant)?;
      slot.node.set_time_distribution(Some(Arc::new(dist_from_parent)));
    }
  }

  Ok(())
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::payload::timetree::NodeTimetree;
  use crate::pretty_assert_ulps_eq;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_distribution::Distribution;

  #[test]
  fn test_forward_pass_set_likely_time_empty_distribution_returns_none() {
    let mut node = node_with_distribution(Some(Distribution::empty()));
    assert_eq!(None, set_likely_time(&mut node, None));
    assert_eq!(None, node.time());
  }

  #[test]
  fn test_forward_pass_set_likely_time_missing_distribution_returns_none() {
    let mut node = node_with_distribution(None);
    assert_eq!(None, set_likely_time(&mut node, None));
    assert_eq!(None, node.time());
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::no_parent(      None,      5.0)]
  #[case::parent_earlier( Some(3.0), 5.0)]
  #[case::parent_clamps(  Some(8.0), 8.0)]
  #[trace]
  fn test_forward_pass_set_likely_time_uses_distribution_peak(
    #[case] parent_time: Option<f64>,
    #[case] expected: f64,
  ) {
    let mut node = node_with_distribution(Some(Distribution::point(5.0, 1.0)));
    let assigned = set_likely_time(&mut node, parent_time).expect("a time should be assigned");
    pretty_assert_ulps_eq!(assigned, expected, max_ulps = 4);
    let committed = node.time().expect("node time should be committed");
    pretty_assert_ulps_eq!(committed, expected, max_ulps = 4);
  }

  mod helpers {
    use super::*;

    pub(super) fn node_with_distribution(distribution: Option<Distribution>) -> NodeTimetree {
      NodeTimetree {
        time_distribution: distribution.map(Arc::new),
        ..NodeTimetree::default()
      }
    }
  }

  use helpers::*;
}
