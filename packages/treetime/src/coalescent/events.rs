use crate::coalescent::time_coordinate::CalendarTime;
use crate::payload::traits::TimetreeNode;
use eyre::Report;
use ordered_float::OrderedFloat;
use treetime_graph::edge::GraphEdge;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNode;
use treetime_utils::make_error;

/// Collects tree merger events as (time, delta_branches) tuples sorted by increasing time.
///
/// Returns present time, events sorted by increasing time (past to present), and
/// the lineage count remaining after the latest retained sample. Each connected
/// bad-branch subtree contributes one remaining lineage because all of its node
/// events are excluded.
/// delta_branches: +1 for leaf nodes, -(k-1) for internal nodes with k children.
pub fn collect_tree_events<N, E, D>(
  graph: &Graph<N, E, D>,
) -> Result<(CalendarTime, Vec<(CalendarTime, i32)>, i32), Report>
where
  N: GraphNode + TimetreeNode,
  E: GraphEdge,
  D: Sync + Send,
{
  if graph.num_roots() != 1 {
    return make_error!("Graph must have exactly one root, found {}", graph.num_roots());
  }

  let mut max_time = CalendarTime::new(f64::NEG_INFINITY);
  let mut events = Vec::new();
  let mut terminal_lineage_count = 0;
  let mut n_excluded = 0;
  let mut n_roots = 0;
  let mut n_leaves = 0;
  let mut n_internal = 0;
  let mut delta_sum = 0;

  graph.iter_breadth_first_forward(|node| {
    let num_children = node.child_edges.len();
    if node.is_leaf != node.child_edges.is_empty() {
      return make_error!(
        "Coalescent event role mismatch for node {:?}: leaf flag is {}, child count is {}",
        node.key,
        node.is_leaf,
        num_children
      );
    }
    if node.is_root != node.parent_keys.is_empty() {
      return make_error!(
        "Coalescent event role mismatch for node {:?}: root flag is {}, parent count is {}",
        node.key,
        node.is_root,
        node.parent_keys.len()
      );
    }

    let delta = if node.is_leaf {
      n_leaves += 1;
      1
    } else {
      n_internal += 1;
      -((num_children as i32) - 1)
    };
    if node.is_root {
      n_roots += 1;
    } else if node.parent_keys.len() != 1 {
      return make_error!(
        "Coalescent events require each non-root node to have exactly one parent, but node {:?} has {}",
        node.key,
        node.parent_keys.len()
      );
    }
    delta_sum += delta;

    if node.payload.bad_branch() {
      terminal_lineage_count += delta;
      n_excluded += 1;
      return Ok(());
    }

    // Every retained node must contribute an event. Silently dropping a good node
    // without an inferred time breaks the lineage-count balance and hides incomplete
    // state after topology changes.
    let Some(t) = node
      .payload
      .time_distribution()
      .as_ref()
      .and_then(|time_dist| time_dist.likely_time())
    else {
      return make_error!(
        "Coalescent lineage count requires an inferred time for every node, but node (key={:?}) has none. \
         The coalescent model was likely built before node times were recomputed for the current tree topology.",
        node.key
      );
    };

    let t = CalendarTime::new(t);
    max_time = max_time.max(t);
    events.push((t, delta));
    Ok(())
  })?;

  if events.is_empty() {
    return make_error!("No tree events found");
  }

  if !max_time.is_finite() {
    return make_error!("Cannot determine present time for coalescent events");
  }

  let expected_nodes = graph.get_nodes().len();
  let expected_leaves = graph.num_leaves();
  let Some(expected_internal) = expected_nodes.checked_sub(expected_leaves) else {
    return make_error!(
      "Coalescent event state is invalid: graph reports {expected_leaves} leaves for {expected_nodes} active nodes"
    );
  };
  if events.len() + n_excluded != expected_nodes
    || n_roots != 1
    || n_leaves != expected_leaves
    || n_internal != expected_internal
    || delta_sum != 1
  {
    return make_error!(
      "Coalescent event state is incomplete: collected {} events for {expected_nodes} nodes, with \
       {n_roots} root(s), {n_leaves}/{expected_leaves} leaves, {n_internal}/{expected_internal} internal nodes, \
       and event delta sum {delta_sum} (expected 1)",
      events.len()
    );
  }

  events.sort_by_key(|(t, _)| OrderedFloat(t.value()));

  Ok((max_time, events, terminal_lineage_count))
}
