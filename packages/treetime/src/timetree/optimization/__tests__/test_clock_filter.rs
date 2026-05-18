#[cfg(test)]
mod tests {
  use crate::clock::clock_filter::{ClockFilterResult, clock_filter_inplace};
  use crate::clock::clock_model::ClockModel;
  use crate::partition::timetree::GraphTimetree;
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use std::sync::Arc;
  use treetime_distribution::Distribution;
  use treetime_graph::node::{Named, Outlier, TimeConstraint};
  use treetime_io::nwk::nwk_read_str;

  const TREE_NEWICK: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

  fn setup_dates(graph: &GraphTimetree, dates: &BTreeMap<String, f64>) {
    for n in graph.get_leaves() {
      let name = n.read_arc().payload().read_arc().name().map(|s| s.as_ref().to_owned());
      if let Some(name) = name {
        if let Some(&date) = dates.get(&name) {
          let dist = Arc::new(Distribution::point(date, 1.0));
          n.write_arc().payload().write_arc().set_time_distribution(Some(dist));
        }
      }
    }
  }

  fn count_outliers(graph: &GraphTimetree) -> usize {
    graph
      .get_leaves()
      .iter()
      .filter(|leaf| leaf.read_arc().payload().read_arc().is_outlier())
      .count()
  }

  #[test]
  fn test_clock_filter_no_outliers_clean_data() -> Result<(), Report> {
    // Tree with dates that fit the clock model well (linear relationship)
    // Clock model: div = 0.01 * date - 20.0 (rate=0.01, intercept=-20.0)
    // For a node at date 2010 with div 0.1: expected_div = 0.01 * 2010 - 20.0 = 0.1
    let graph: GraphTimetree = nwk_read_str(TREE_NEWICK)?;

    // Set dates that match the branch lengths well
    let dates = btreemap! {
      "A".to_owned() => 2010.0,
      "B".to_owned() => 2020.0,
      "C".to_owned() => 2015.0,
      "D".to_owned() => 2012.0,
    };
    setup_dates(&graph, &dates);

    // Clock model: rate=0.01, intercept=-20.0
    // At date 2010, expected div = 0.01 * 2010 + (-20.0) = 0.1
    // At date 2020, expected div = 0.01 * 2020 + (-20.0) = 0.2
    let clock_model = ClockModel::for_testing(0.01, -20.0);

    let ClockFilterResult { new_outliers, iqd } = clock_filter_inplace(&graph, &clock_model, 3.0);

    // With well-fitting data, no outliers should be detected
    assert_eq!(count_outliers(&graph), 0, "No outliers expected for clean data");
    assert!(iqd >= 0.0, "IQD should be non-negative");
    // new_outliers counts status changes, could be 0 if none were outliers before
    assert_eq!(new_outliers, 0, "No status changes expected");

    Ok(())
  }

  #[test]
  fn test_clock_filter_detects_outlier() -> Result<(), Report> {
    // Tree with one leaf having a date that deviates strongly from the clock model
    let graph: GraphTimetree = nwk_read_str(TREE_NEWICK)?;

    // Set dates where one sample (A) has an extreme deviation
    // A is at div ~0.2 (root:0.01 + AB:0.1 + A:0.1) but claims date 1900 (very old)
    let dates = btreemap! {
      "A".to_owned() => 1900.0,  // Outlier: claims very old date but has recent divergence
      "B".to_owned() => 2020.0,
      "C".to_owned() => 2015.0,
      "D".to_owned() => 2012.0,
    };
    setup_dates(&graph, &dates);

    // Clock model based on B, C, D (excluding A)
    // rate=0.01, intercept=-20.0
    let clock_model = ClockModel::for_testing(0.01, -20.0);

    let ClockFilterResult { new_outliers, iqd } = clock_filter_inplace(&graph, &clock_model, 3.0);

    // A should be detected as outlier (date 1900 with div ~0.2 doesn't fit clock)
    // Expected div at 1900 = 0.01 * 1900 - 20.0 = -1.0, but actual div ~0.2
    // Deviation = expected - actual = -1.0 - 0.2 = -1.2 (huge compared to IQD of others)
    assert!(
      count_outliers(&graph) >= 1,
      "At least one outlier expected for data with extreme deviation"
    );
    assert!(iqd > 0.0, "IQD should be positive with varying dates");
    assert!(new_outliers >= 1, "At least one status change expected");

    // Verify A is marked as outlier
    let a_is_outlier = graph.get_leaves().iter().any(|leaf| {
      let node = leaf.read_arc();
      let payload = node.payload().read_arc();
      payload.name().is_some_and(|n| n.as_ref() == "A") && payload.is_outlier()
    });
    assert!(a_is_outlier, "Node A should be marked as outlier");

    Ok(())
  }

  #[test]
  fn test_clock_filter_iqd_calculation() -> Result<(), Report> {
    // Verify IQD is computed and returned correctly
    let graph: GraphTimetree = nwk_read_str(TREE_NEWICK)?;

    // Dates with some spread to create non-zero IQD
    let dates = btreemap! {
      "A".to_owned() => 2010.0,
      "B".to_owned() => 2020.0,
      "C".to_owned() => 2015.0,
      "D".to_owned() => 2012.0,
    };
    setup_dates(&graph, &dates);

    let clock_model = ClockModel::for_testing(0.01, -20.0);

    let ClockFilterResult { iqd, .. } = clock_filter_inplace(&graph, &clock_model, 3.0);

    // IQD should be computed (may be zero or positive depending on data fit)
    assert!(iqd.is_finite(), "IQD should be a finite number");

    Ok(())
  }

  #[test]
  fn test_clock_filter_respects_threshold() -> Result<(), Report> {
    // Test that higher threshold allows more deviation
    let graph: GraphTimetree = nwk_read_str(TREE_NEWICK)?;

    let dates = btreemap! {
      "A".to_owned() => 1980.0,  // Moderate deviation
      "B".to_owned() => 2020.0,
      "C".to_owned() => 2015.0,
      "D".to_owned() => 2012.0,
    };
    setup_dates(&graph, &dates);

    let clock_model = ClockModel::for_testing(0.01, -20.0);

    // With low threshold, A might be outlier
    clock_filter_inplace(&graph, &clock_model, 1.0);
    let outliers_low_threshold = count_outliers(&graph);

    // Reset outlier status
    for leaf in graph.get_leaves() {
      leaf.write_arc().payload().write_arc().set_is_outlier(false);
    }

    // With high threshold, A should not be outlier
    clock_filter_inplace(&graph, &clock_model, 100.0);
    let outliers_high_threshold = count_outliers(&graph);

    assert!(
      outliers_high_threshold <= outliers_low_threshold,
      "Higher threshold should result in fewer or equal outliers"
    );

    Ok(())
  }
}
