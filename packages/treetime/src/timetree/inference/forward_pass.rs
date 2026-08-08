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
  refine_distribution_from_parent(dependencies, slot)?;

  let has_exact_date = has_exact_date(&slot.node);

  // Independent marginal modes can invert even though the forward pass commits each parent before
  // visiting its children. Current behavior projects the point estimate of a node whose time was
  // inferred to the committed parent time; a node given an exact date keeps it, so that a date
  // conflicting with the fitted clock stays visible to `commit_clock_branch_lengths` rather than
  // being clamped away. This changes the point estimate without recomputing the posterior. The
  // statistical contract remains open in
  // kb/issues/M-timetree-marginal-node-times-can-violate-topology.md.
  let parent_time = (!has_exact_date).then(|| parent_time(dependencies, slot)).flatten();

  // A node whose time distribution is empty (or degenerate) yields no likely time, so no date is
  // assigned. This signals irreconcilable inference constraints at the node (e.g. disjoint child
  // messages, or a date range disjoint from what the rest of the tree implies). Surface it instead
  // of dropping the date silently. An undated leaf under an undated parent has nothing to infer
  // from and is not worth reporting: it was already accounted for when its date was found missing.
  let is_dateable = !graph.is_leaf(slot.key) || slot.node.date_constraint().is_some();
  if set_likely_time(&mut slot.node, parent_time).is_none() && is_dateable {
    let name = slot.node.name();
    let name = name.as_ref().map_or("<unnamed>", |name| name.as_ref());
    warn!(
      "Timetree forward pass: node '{name}' has an empty time distribution; no date was \
       assigned. The inference constraints at this node are likely irreconcilable (e.g. conflicting \
       or disjoint date constraints from its subtree)."
    );
  }
  Ok(())
}

/// Whether the node was given an exact date, which pins it to a single time.
///
/// Anything else -- an uncertain or ranged date, or no date at all -- leaves the node's time to be
/// inferred, and so to be refined by the message coming down from its parent.
fn has_exact_date(node: &impl TimetreeNode) -> bool {
  node.date_constraint().as_ref().is_some_and(|dist| dist.is_point())
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

fn refine_distribution_from_parent<N, E>(
  dependencies: &IndexedPassDependencies<N, E>,
  slot: &mut IndexedPassSlot<N, E>,
) -> Result<(), Report>
where
  N: GraphNode + TimetreeNode,
  E: GraphEdge + TimetreeEdge,
{
  if slot.parent_key.is_none() {
    return Ok(());
  }

  // An exactly dated node has nothing to refine: its time is observed, not inferred. Leaves whose
  // date is uncertain or missing are refined like any other node -- the message from the parent is
  // the only thing that can narrow an ambiguous date down.
  if has_exact_date(&slot.node) {
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
      // Same forward tail policy as the two-parent-message branch above. The forward message
      // needs Zero left / Constant right regardless of the convolution's own tails, so apply
      // them explicitly after normalize() (which preserves whatever tails it received).
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
