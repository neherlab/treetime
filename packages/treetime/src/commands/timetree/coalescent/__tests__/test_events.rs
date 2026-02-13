#[cfg(test)]
mod tests {
  use crate::commands::clock::date_constraints::load_date_constraints;
  use crate::commands::timetree::coalescent::events::collect_tree_events;
  use crate::representation::partition_timetree::GraphTimetree;
  use approx::assert_ulps_eq;
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use treetime_io::dates_csv::{DateOrRange, DatesMap};
  use treetime_io::nwk::nwk_read_str;

  fn create_graph_with_dates(tree_nwk: &str, dates: &DatesMap) -> Result<GraphTimetree, Report> {
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

    assert_ulps_eq!(present_time, 2015.0, max_ulps = 4);
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

    assert_ulps_eq!(present_time, 2012.0, max_ulps = 4);
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

    assert_ulps_eq!(present_time, 2010.0, max_ulps = 4);
    assert_eq!(events, vec![(2000.0, -2), (2010.0, 1), (2010.0, 1), (2010.0, 1)]);

    Ok(())
  }
}
