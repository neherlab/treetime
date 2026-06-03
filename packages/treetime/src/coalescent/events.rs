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
/// Returns present time and events sorted by increasing time (past to present).
/// delta_branches: +1 for leaf nodes, -(k-1) for internal nodes with k children.
pub fn collect_tree_events<N, E, D>(graph: &Graph<N, E, D>) -> Result<(CalendarTime, Vec<(CalendarTime, i32)>), Report>
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

  graph.iter_breadth_first_forward(|node| {
    if let Some(time_dist) = node.payload.time_distribution() {
      if let Some(t) = time_dist.likely_time() {
        let t = CalendarTime::new(t);
        max_time = max_time.max(t);

        let num_children = node.child_edges.len();

        if num_children == 0 {
          events.push((t, 1));
        } else {
          events.push((t, -((num_children as i32) - 1)));
        }
      }
    }
    Ok(())
  })?;

  if events.is_empty() {
    return make_error!("No tree events found");
  }

  if !max_time.is_finite() {
    return make_error!("Cannot determine present time for coalescent events");
  }

  events.sort_by_key(|(t, _)| OrderedFloat(t.value()));

  Ok((max_time, events))
}
