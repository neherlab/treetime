#[cfg(test)]
mod tests {
  use crate::clock::date_constraints::load_date_constraints;
  use crate::coalescent::events::collect_tree_events;
  use crate::coalescent::time_coordinate::CalendarTime;
  use crate::partition::timetree::partition::GraphTimetree;
  use crate::payload::timetree::NodeTimetree;
  use crate::pretty_assert_ulps_eq;
  use crate::test_utils::find_node_key_by_name;
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use std::sync::Arc;
  use treetime_distribution::Distribution;
  use treetime_io::dates_csv::{DateConstraint, DatesMap};
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::assert_error;

  fn create_graph_with_dates(tree_nwk: &str, dates: &DatesMap) -> Result<GraphTimetree, Report> {
    let graph = nwk_read_str(tree_nwk)?;
    load_date_constraints(dates, &graph)?;
    Ok(graph)
  }

  fn cal(t: f64) -> CalendarTime {
    CalendarTime::new(t)
  }

  #[test]
  fn test_collect_tree_events_simple() -> Result<(), Report> {
    const TREE_NWK: &str = "(child1:1.0,child2:1.0,child3:1.0)root:1.0;";
    let dates = btreemap! {
      "root".to_owned() => Some(DateConstraint::exact(2000.0)),
      "child1".to_owned() => Some(DateConstraint::exact(2005.0)),
      "child2".to_owned() => Some(DateConstraint::exact(2010.0)),
      "child3".to_owned() => Some(DateConstraint::exact(2015.0)),
    };

    let graph = create_graph_with_dates(TREE_NWK, &dates)?;
    let (present_time, events, terminal_lineage_count) = collect_tree_events(&graph)?;

    pretty_assert_ulps_eq!(present_time.value(), 2015.0, max_ulps = 4);
    assert_eq!(terminal_lineage_count, 0);
    assert_eq!(
      events,
      vec![(cal(2000.0), -2), (cal(2005.0), 1), (cal(2010.0), 1), (cal(2015.0), 1)]
    );

    Ok(())
  }

  #[test]
  fn test_collect_tree_events_complex() -> Result<(), Report> {
    const TREE_NWK: &str = "((leaf1:1.0,leaf2:1.0)internal1:1.0,leaf3:1.0)root:1.0;";
    let dates = btreemap! {
      "root".to_owned() => Some(DateConstraint::exact(2000.0)),
      "internal1".to_owned() => Some(DateConstraint::exact(2005.0)),
      "leaf1".to_owned() => Some(DateConstraint::exact(2010.0)),
      "leaf2".to_owned() => Some(DateConstraint::exact(2010.0)),
      "leaf3".to_owned() => Some(DateConstraint::exact(2012.0)),
    };

    let graph = create_graph_with_dates(TREE_NWK, &dates)?;
    let (present_time, events, terminal_lineage_count) = collect_tree_events(&graph)?;

    pretty_assert_ulps_eq!(present_time.value(), 2012.0, max_ulps = 4);
    assert_eq!(terminal_lineage_count, 0);
    assert_eq!(
      events,
      vec![
        (cal(2000.0), -1),
        (cal(2005.0), -1),
        (cal(2010.0), 1),
        (cal(2010.0), 1),
        (cal(2012.0), 1)
      ]
    );

    Ok(())
  }

  #[test]
  fn test_collect_tree_events_sorted() -> Result<(), Report> {
    const TREE_NWK: &str = "((leaf1:1.0,leaf2:1.0)internal1:1.0,leaf3:1.0)root:1.0;";
    let dates = btreemap! {
      "root".to_owned() => Some(DateConstraint::exact(2000.0)),
      "internal1".to_owned() => Some(DateConstraint::exact(2005.0)),
      "leaf1".to_owned() => Some(DateConstraint::exact(2010.0)),
      "leaf2".to_owned() => Some(DateConstraint::exact(2010.0)),
      "leaf3".to_owned() => Some(DateConstraint::exact(2012.0)),
    };

    let graph = create_graph_with_dates(TREE_NWK, &dates)?;
    let (_present_time, events, _terminal_lineage_count) = collect_tree_events(&graph)?;

    for i in 1..events.len() {
      assert!(events[i - 1].0 <= events[i].0);
    }

    Ok(())
  }

  #[test]
  fn test_collect_tree_events_rejects_node_without_time() -> Result<(), Report> {
    // internal1 has no date constraint, so its time distribution stays unset and it
    // cannot contribute its merger event. Event collection must reject this rather
    // than silently drop the node (which would unbalance the lineage-count deltas).
    const TREE_NWK: &str = "((leaf1:1.0,leaf2:1.0)internal1:1.0,leaf3:1.0)root:1.0;";
    let dates = btreemap! {
      "root".to_owned() => Some(DateConstraint::exact(2000.0)),
      "leaf1".to_owned() => Some(DateConstraint::exact(2010.0)),
      "leaf2".to_owned() => Some(DateConstraint::exact(2010.0)),
      "leaf3".to_owned() => Some(DateConstraint::exact(2012.0)),
    };

    let graph = create_graph_with_dates(TREE_NWK, &dates)?;

    assert_error!(
      collect_tree_events(&graph),
      "Coalescent lineage count requires an inferred time for every node, but node (key=GraphNodeKey(2)) has none. The coalescent model was likely built before node times were recomputed for the current tree topology."
    );

    Ok(())
  }

  #[test]
  fn test_collect_tree_events_rejects_unreachable_active_node() -> Result<(), Report> {
    const TREE_NWK: &str = "(child1:1.0,child2:1.0,child3:1.0)root:1.0;";
    let dates = btreemap! {
      "root".to_owned() => Some(DateConstraint::exact(2000.0)),
      "child1".to_owned() => Some(DateConstraint::exact(2005.0)),
      "child2".to_owned() => Some(DateConstraint::exact(2010.0)),
      "child3".to_owned() => Some(DateConstraint::exact(2015.0)),
    };
    let mut graph = create_graph_with_dates(TREE_NWK, &dates)?;
    graph.add_node(NodeTimetree {
      time: Some(2002.0),
      time_distribution: Some(Arc::new(Distribution::point(2002.0, 1.0))),
      ..NodeTimetree::default()
    });

    assert_error!(
      collect_tree_events(&graph),
      "Coalescent event state is incomplete: collected 4 events for 5 nodes, with 1 root(s), 3/3 leaves, 1/2 internal nodes, and event delta sum 1 (expected 1)"
    );

    Ok(())
  }

  #[test]
  fn test_collect_tree_events_simultaneous_events() -> Result<(), Report> {
    const TREE_NWK: &str = "(child1:1.0,child2:1.0,child3:1.0)root:1.0;";
    let dates = btreemap! {
      "root".to_owned() => Some(DateConstraint::exact(2000.0)),
      "child1".to_owned() => Some(DateConstraint::exact(2010.0)),
      "child2".to_owned() => Some(DateConstraint::exact(2010.0)),
      "child3".to_owned() => Some(DateConstraint::exact(2010.0)),
    };

    let graph = create_graph_with_dates(TREE_NWK, &dates)?;
    let (present_time, events, terminal_lineage_count) = collect_tree_events(&graph)?;

    pretty_assert_ulps_eq!(present_time.value(), 2010.0, max_ulps = 4);
    assert_eq!(terminal_lineage_count, 0);
    assert_eq!(
      events,
      vec![(cal(2000.0), -2), (cal(2010.0), 1), (cal(2010.0), 1), (cal(2010.0), 1)]
    );

    Ok(())
  }

  #[test]
  fn test_collect_tree_events_excludes_bad_leaf() -> Result<(), Report> {
    const TREE_NWK: &str = "(child1:1.0,child2:1.0,child3:1.0)root:1.0;";
    let dates = btreemap! {
      "root".to_owned() => Some(DateConstraint::exact(2000.0)),
      "child1".to_owned() => Some(DateConstraint::exact(2005.0)),
      "child2".to_owned() => Some(DateConstraint::exact(2010.0)),
      "child3".to_owned() => Some(DateConstraint::exact(2015.0)),
    };
    let graph = create_graph_with_dates(TREE_NWK, &dates)?;
    let child3_key = find_node_key_by_name(&graph, "child3").expect("child3 not found");
    graph
      .get_node(child3_key)
      .expect("child3 exists")
      .read_arc()
      .payload()
      .write_arc()
      .bad_branch = true;

    let (present_time, events, terminal_lineage_count) = collect_tree_events(&graph)?;

    // Oracle: v0 filters `bad_branch` nodes before constructing `tree_events`.
    // packages/legacy/treetime/treetime/merger_models.py#L102-L105
    pretty_assert_ulps_eq!(present_time.value(), 2010.0, max_ulps = 4);
    assert_eq!(terminal_lineage_count, 1);
    assert_eq!(events, vec![(cal(2000.0), -2), (cal(2005.0), 1), (cal(2010.0), 1)]);

    Ok(())
  }

  #[test]
  fn test_collect_tree_events_excludes_bad_leaf_without_time() -> Result<(), Report> {
    const TREE_NWK: &str = "(child1:1.0,child2:1.0,child3:1.0,child4:1.0)root:1.0;";
    let dates = btreemap! {
      "root".to_owned() => Some(DateConstraint::exact(2000.0)),
      "child1".to_owned() => Some(DateConstraint::exact(2005.0)),
      "child2".to_owned() => Some(DateConstraint::exact(2010.0)),
      "child3".to_owned() => Some(DateConstraint::exact(2015.0)),
    };
    let graph = create_graph_with_dates(TREE_NWK, &dates)?;

    let (present_time, events, terminal_lineage_count) = collect_tree_events(&graph)?;

    pretty_assert_ulps_eq!(present_time.value(), 2015.0, max_ulps = 4);
    assert_eq!(terminal_lineage_count, 1);
    assert_eq!(
      events,
      vec![(cal(2000.0), -3), (cal(2005.0), 1), (cal(2010.0), 1), (cal(2015.0), 1)]
    );

    Ok(())
  }

  #[test]
  fn test_collect_tree_events_counts_excluded_subtree_once() -> Result<(), Report> {
    const TREE_NWK: &str = "((leaf1:1.0,leaf2:1.0)internal1:1.0,leaf3:1.0,leaf4:1.0,leaf5:1.0)root:1.0;";
    let dates = btreemap! {
      "root".to_owned() => Some(DateConstraint::exact(2000.0)),
      "leaf3".to_owned() => Some(DateConstraint::exact(2005.0)),
      "leaf4".to_owned() => Some(DateConstraint::exact(2010.0)),
      "leaf5".to_owned() => Some(DateConstraint::exact(2015.0)),
    };
    let graph = create_graph_with_dates(TREE_NWK, &dates)?;

    let (present_time, events, terminal_lineage_count) = collect_tree_events(&graph)?;

    pretty_assert_ulps_eq!(present_time.value(), 2015.0, max_ulps = 4);
    assert_eq!(terminal_lineage_count, 1);
    assert_eq!(
      events,
      vec![(cal(2000.0), -3), (cal(2005.0), 1), (cal(2010.0), 1), (cal(2015.0), 1)]
    );

    Ok(())
  }
}
