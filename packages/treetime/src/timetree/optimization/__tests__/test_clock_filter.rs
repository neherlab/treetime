#[cfg(test)]
mod tests {
  use crate::clock::clock_filter::{ClockFilterResult, clock_filter_inplace};
  use crate::clock::clock_model::ClockModel;
  use crate::clock::clock_state::ClockState;
  use crate::clock::date_constraints::DateConstraints;
  use crate::partition::timetree::partition::GraphTimetree;
  use crate::timetree::optimization::clock_filter::propagate_bad_branches;
  use crate::timetree::timetree_state::TimetreeState;
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use std::sync::Arc;
  use treetime_distribution::Distribution;
  use treetime_graph::node::Named;
  use treetime_io::nwk::nwk_read_str;

  const TREE_NEWICK: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

  /// Per-leaf date inputs as a value map, replacing the payload `time_distribution` the tests used to
  /// write. Each named leaf with a date carries a point distribution at that date.
  fn date_constraints(graph: &GraphTimetree, dates: &BTreeMap<String, f64>) -> DateConstraints {
    let mut time_distributions = BTreeMap::new();
    for n in graph.get_leaves() {
      let n = n.read_arc();
      let name = n.payload().read_arc().name().map(|s| s.as_ref().to_owned());
      if let Some(name) = name {
        if let Some(&date) = dates.get(&name) {
          time_distributions.insert(n.key(), Some(Arc::new(Distribution::point(date, 1.0))));
        }
      }
    }
    DateConstraints {
      time_distributions,
      ..DateConstraints::default()
    }
  }

  /// Seed the clock state from date-constraint values, reproducing the payload-reading seed: the date
  /// state built from `constraints` supplies each node's date through `likely_times`, and the clock
  /// state starts from those dates with default divergence and outlier flags.
  fn seed_clock_state(graph: &GraphTimetree, constraints: &DateConstraints) -> ClockState {
    let date_state = TimetreeState::seed_from_values(graph, constraints);
    ClockState::seed_from_values(graph, &date_state.likely_times())
  }

  fn count_outliers(graph: &GraphTimetree, state: &ClockState) -> usize {
    graph
      .get_leaves()
      .iter()
      .filter(|leaf| state.node(leaf.read_arc().key()).is_outlier)
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
    let constraints = date_constraints(&graph, &dates);

    // Clock model: rate=0.01, intercept=-20.0
    // At date 2010, expected div = 0.01 * 2010 + (-20.0) = 0.1
    // At date 2020, expected div = 0.01 * 2020 + (-20.0) = 0.2
    let clock_model = ClockModel::for_testing(0.01, -20.0);

    let mut state = seed_clock_state(&graph, &constraints);
    let ClockFilterResult { new_outliers, iqd } = clock_filter_inplace(&graph, &mut state, &clock_model, 3.0)?;

    // With well-fitting data, no outliers should be detected
    assert_eq!(count_outliers(&graph, &state), 0, "No outliers expected for clean data");
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
    let constraints = date_constraints(&graph, &dates);

    // Clock model based on B, C, D (excluding A)
    // rate=0.01, intercept=-20.0
    let clock_model = ClockModel::for_testing(0.01, -20.0);

    let mut state = seed_clock_state(&graph, &constraints);
    let ClockFilterResult { new_outliers, iqd } = clock_filter_inplace(&graph, &mut state, &clock_model, 3.0)?;

    // A should be detected as outlier (date 1900 with div ~0.2 doesn't fit clock)
    // Expected div at 1900 = 0.01 * 1900 - 20.0 = -1.0, but actual div ~0.2
    // Deviation = expected - actual = -1.0 - 0.2 = -1.2 (huge compared to IQD of others)
    assert!(
      count_outliers(&graph, &state) >= 1,
      "At least one outlier expected for data with extreme deviation"
    );
    assert!(iqd > 0.0, "IQD should be positive with varying dates");
    assert!(new_outliers >= 1, "At least one status change expected");

    // Verify A is marked as outlier
    let a_is_outlier = graph.get_leaves().iter().any(|leaf| {
      let node = leaf.read_arc();
      state.node(node.key()).is_outlier && node.payload().read_arc().name().is_some_and(|n| n.as_ref() == "A")
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
    let constraints = date_constraints(&graph, &dates);

    let clock_model = ClockModel::for_testing(0.01, -20.0);

    let mut state = seed_clock_state(&graph, &constraints);
    let ClockFilterResult { iqd, .. } = clock_filter_inplace(&graph, &mut state, &clock_model, 3.0)?;

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
    let constraints = date_constraints(&graph, &dates);

    let clock_model = ClockModel::for_testing(0.01, -20.0);

    // With low threshold, A might be outlier
    let mut state_low = seed_clock_state(&graph, &constraints);
    clock_filter_inplace(&graph, &mut state_low, &clock_model, 1.0)?;
    let outliers_low_threshold = count_outliers(&graph, &state_low);

    // With high threshold, A should not be outlier. Each filter runs on its own freshly seeded state,
    // so the low-threshold outlier flags do not carry over.
    let mut state_high = seed_clock_state(&graph, &constraints);
    clock_filter_inplace(&graph, &mut state_high, &clock_model, 100.0)?;
    let outliers_high_threshold = count_outliers(&graph, &state_high);

    assert!(
      outliers_high_threshold <= outliers_low_threshold,
      "Higher threshold should result in fewer or equal outliers"
    );

    Ok(())
  }

  #[test]
  fn test_clock_filter_propagates_bad_branches_after_topology_change() -> Result<(), Report> {
    let graph: GraphTimetree = nwk_read_str("((A:0.1,B:0.1)AB:0.1,C:0.1)root;")?;
    let mut state = TimetreeState::new(&graph);
    for node in graph.get_leaves() {
      let node = node.read_arc();
      let is_bad = node
        .payload()
        .read_arc()
        .name()
        .is_some_and(|name| name.as_ref() != "C");
      state.node_mut(node.key()).bad_branch = is_bad;
    }

    propagate_bad_branches(&graph, &mut state)?;

    let actual = graph
      .get_nodes()
      .iter()
      .map(|node| {
        let node = node.read_arc();
        let name = node
          .payload()
          .read_arc()
          .name()
          .expect("Every fixture node must be named")
          .as_ref()
          .to_owned();
        (name, state.node(node.key()).bad_branch)
      })
      .collect::<BTreeMap<_, _>>();
    let expected = btreemap! {
      "A".to_owned() => true,
      "AB".to_owned() => true,
      "B".to_owned() => true,
      "C".to_owned() => false,
      "root".to_owned() => false,
    };
    assert_eq!(expected, actual);

    Ok(())
  }
}
