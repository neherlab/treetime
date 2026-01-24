use crate::commands::timetree::timetree_traits::TimetreeNode;
use crate::graph::edge::GraphEdge;
use crate::graph::graph::Graph;
use crate::graph::node::GraphNode;
use eyre::Report;
use ordered_float::OrderedFloat;
use treetime_utils::make_error;

/// Collects tree merger events as (time, delta_branches) tuples sorted by increasing time.
///
/// Returns present time and events sorted by increasing time (past to present).
/// delta_branches: +1 for leaf nodes, -(k-1) for internal nodes with k children.
pub fn collect_tree_events<N, E, D>(graph: &Graph<N, E, D>) -> Result<(f64, Vec<(f64, i32)>), Report>
where
  N: GraphNode + TimetreeNode,
  E: GraphEdge,
  D: Sync + Send,
{
  if graph.num_roots() != 1 {
    return make_error!("Graph must have exactly one root, found {}", graph.num_roots());
  }

  let mut max_time = f64::NEG_INFINITY;
  let mut events = Vec::new();

  graph.iter_breadth_first_forward(|node| {
    if let Some(time_dist) = node.payload.time_distribution() {
      if let Some(t) = time_dist.likely_time() {
        max_time = max_time.max(t);

        let num_children = node.child_edges.len();

        if num_children == 0 {
          events.push((t, 1));
        } else {
          events.push((t, -((num_children as i32) - 1)));
        }
      }
    }
  });

  if events.is_empty() {
    return make_error!("No tree events found");
  }

  if !max_time.is_finite() {
    return make_error!("Cannot determine present time for coalescent events");
  }

  events.sort_by_key(|x| OrderedFloat(x.0));

  Ok((max_time, events))
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::commands::timetree::data::date_constraints::load_date_constraints;
  use crate::io::dates_csv::{DateOrRange, DatesMap};
  use crate::io::nwk::nwk_read_str;
  use crate::representation::graph_ancestral::GraphAncestral;
  use approx::assert_abs_diff_eq;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;

  fn create_graph_with_dates(tree_nwk: &str, dates: &DatesMap) -> Result<GraphAncestral, Report> {
    let graph = nwk_read_str(tree_nwk)?;
    load_date_constraints(dates, &graph)?;
    Ok(graph)
  }

  #[test]
  fn test_collect_tree_events_simple() -> Result<(), Report> {
    const TREE_NWK: &str = "(child1:1.0,child2:1.0,child3:1.0)root:1.0;";
    let dates = btreemap! {
      "root".to_owned() => Some(DateOrRange::YearFraction(2000.0)),
      "child1".to_owned() => Some(DateOrRange::YearFraction(2005.0)),
      "child2".to_owned() => Some(DateOrRange::YearFraction(2010.0)),
      "child3".to_owned() => Some(DateOrRange::YearFraction(2015.0)),
    };

    let graph = create_graph_with_dates(TREE_NWK, &dates)?;
    let (present_time, events) = collect_tree_events(&graph)?;

    assert_abs_diff_eq!(present_time, 2015.0, epsilon = 1e-10);
    assert_eq!(events, vec![(2000.0, -2), (2005.0, 1), (2010.0, 1), (2015.0, 1)]);

    Ok(())
  }

  #[test]
  fn test_collect_tree_events_complex() -> Result<(), Report> {
    const TREE_NWK: &str = "((leaf1:1.0,leaf2:1.0)internal1:1.0,leaf3:1.0)root:1.0;";
    let dates = btreemap! {
      "root".to_owned() => Some(DateOrRange::YearFraction(2000.0)),
      "internal1".to_owned() => Some(DateOrRange::YearFraction(2005.0)),
      "leaf1".to_owned() => Some(DateOrRange::YearFraction(2010.0)),
      "leaf2".to_owned() => Some(DateOrRange::YearFraction(2010.0)),
      "leaf3".to_owned() => Some(DateOrRange::YearFraction(2012.0)),
    };

    let graph = create_graph_with_dates(TREE_NWK, &dates)?;
    let (present_time, events) = collect_tree_events(&graph)?;

    assert_abs_diff_eq!(present_time, 2012.0, epsilon = 1e-10);
    assert_eq!(
      events,
      vec![(2000.0, -1), (2005.0, -1), (2010.0, 1), (2010.0, 1), (2012.0, 1)]
    );

    Ok(())
  }

  #[test]
  fn test_collect_tree_events_sorted() -> Result<(), Report> {
    const TREE_NWK: &str = "((leaf1:1.0,leaf2:1.0)internal1:1.0,leaf3:1.0)root:1.0;";
    let dates = btreemap! {
      "root".to_owned() => Some(DateOrRange::YearFraction(2000.0)),
      "internal1".to_owned() => Some(DateOrRange::YearFraction(2005.0)),
      "leaf1".to_owned() => Some(DateOrRange::YearFraction(2010.0)),
      "leaf2".to_owned() => Some(DateOrRange::YearFraction(2010.0)),
      "leaf3".to_owned() => Some(DateOrRange::YearFraction(2012.0)),
    };

    let graph = create_graph_with_dates(TREE_NWK, &dates)?;
    let (_present_time, events) = collect_tree_events(&graph)?;

    for i in 1..events.len() {
      assert!(events[i - 1].0 <= events[i].0);
    }

    Ok(())
  }

  #[test]
  fn test_collect_tree_events_simultaneous_events() -> Result<(), Report> {
    const TREE_NWK: &str = "(child1:1.0,child2:1.0,child3:1.0)root:1.0;";
    let dates = btreemap! {
      "root".to_owned() => Some(DateOrRange::YearFraction(2000.0)),
      "child1".to_owned() => Some(DateOrRange::YearFraction(2010.0)),
      "child2".to_owned() => Some(DateOrRange::YearFraction(2010.0)),
      "child3".to_owned() => Some(DateOrRange::YearFraction(2010.0)),
    };

    let graph = create_graph_with_dates(TREE_NWK, &dates)?;
    let (present_time, events) = collect_tree_events(&graph)?;

    assert_abs_diff_eq!(present_time, 2010.0, epsilon = 1e-10);
    assert_eq!(events, vec![(2000.0, -2), (2010.0, 1), (2010.0, 1), (2010.0, 1)]);

    Ok(())
  }
}
