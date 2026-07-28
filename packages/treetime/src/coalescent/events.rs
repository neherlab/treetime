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

  graph.iter_breadth_first_forward(|node| {
    let delta = tree_event_delta(node.child_edges.len());
    if node.payload.bad_branch() {
      terminal_lineage_count += delta;
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

  events.sort_by_key(|(t, _)| OrderedFloat(t.value()));

  Ok((max_time, events, terminal_lineage_count))
}

fn tree_event_delta(num_children: usize) -> i32 {
  if num_children == 0 {
    1
  } else {
    -((num_children as i32) - 1)
  }
}
