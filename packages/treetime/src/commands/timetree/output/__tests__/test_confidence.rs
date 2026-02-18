#[cfg(test)]
mod tests {
  use crate::commands::timetree::output::confidence::{combine_confidence, extract_confidence_intervals};
  use crate::representation::partition::timetree::GraphTimetree;
  use crate::representation::payload::timetree::NodeTimetree;
  use approx::assert_relative_eq;
  use std::sync::Arc;
  use treetime_distribution::Distribution;
  use treetime_graph::node::Named;

  fn make_node(name: Option<&str>, time: Option<f64>, dist: Option<Arc<Distribution>>) -> NodeTimetree {
    let mut node = NodeTimetree::default();
    node.base.set_name(name);
    node.time = time;
    node.time_distribution = dist;
    node
  }

  #[test]
  fn test_extract_confidence_intervals_skips_unnamed_nodes() {
    let mut graph = GraphTimetree::new();
    graph.add_node(make_node(None, Some(2020.0), None));
    graph.add_node(make_node(Some("named"), Some(2021.0), None));
    graph.build().unwrap();

    let intervals = extract_confidence_intervals(&graph);
    assert_eq!(intervals.len(), 1);
    assert_eq!(intervals[0].name, "named");
  }

  #[test]
  fn test_extract_confidence_intervals_skips_nodes_without_time() {
    let mut graph = GraphTimetree::new();
    graph.add_node(make_node(Some("no_time"), None, None));
    graph.add_node(make_node(Some("has_time"), Some(2021.0), None));
    graph.build().unwrap();

    let intervals = extract_confidence_intervals(&graph);
    assert_eq!(intervals.len(), 1);
    assert_eq!(intervals[0].name, "has_time");
  }

  #[test]
  fn test_extract_confidence_intervals_uses_date_as_fallback() {
    let mut graph = GraphTimetree::new();
    graph.add_node(make_node(Some("node_a"), Some(2020.5), None));
    graph.build().unwrap();

    let intervals = extract_confidence_intervals(&graph);
    assert_eq!(intervals.len(), 1);
    assert_relative_eq!(intervals[0].date, 2020.5);
    assert_relative_eq!(intervals[0].lower, 2020.5);
    assert_relative_eq!(intervals[0].upper, 2020.5);
  }

  #[test]
  fn test_extract_confidence_intervals_with_distribution() {
    let mut graph = GraphTimetree::new();
    let dist = Arc::new(Distribution::range((2019.0, 2021.0), 1.0));
    graph.add_node(make_node(Some("node_a"), Some(2020.0), Some(dist)));
    graph.build().unwrap();

    let intervals = extract_confidence_intervals(&graph);
    assert_eq!(intervals.len(), 1);
    assert_relative_eq!(intervals[0].date, 2020.0);
    // 95% CI from uniform: 0.025 * 2 + 2019 = 2019.05, 0.975 * 2 + 2019 = 2020.95
    assert_relative_eq!(intervals[0].lower, 2019.05, epsilon = 1e-10);
    assert_relative_eq!(intervals[0].upper, 2020.95, epsilon = 1e-10);
  }

  #[test]
  fn test_extract_confidence_intervals_sorted_by_name() {
    let mut graph = GraphTimetree::new();
    graph.add_node(make_node(Some("zebra"), Some(2020.0), None));
    graph.add_node(make_node(Some("alpha"), Some(2021.0), None));
    graph.add_node(make_node(Some("middle"), Some(2022.0), None));
    graph.build().unwrap();

    let intervals = extract_confidence_intervals(&graph);
    assert_eq!(intervals.len(), 3);
    assert_eq!(intervals[0].name, "alpha");
    assert_eq!(intervals[1].name, "middle");
    assert_eq!(intervals[2].name, "zebra");
  }

  #[test]
  fn test_combine_confidence_no_contributions() {
    let result = combine_confidence(10.0, (5.0, 15.0), None, None);
    assert_relative_eq!(result.0, 5.0);
    assert_relative_eq!(result.1, 15.0);
  }

  #[test]
  fn test_combine_confidence_single_contribution() {
    let result = combine_confidence(10.0, (0.0, 20.0), Some((8.0, 12.0)), None);
    assert_relative_eq!(result.0, 8.0);
    assert_relative_eq!(result.1, 12.0);

    let result = combine_confidence(10.0, (0.0, 20.0), None, Some((7.0, 13.0)));
    assert_relative_eq!(result.0, 7.0);
    assert_relative_eq!(result.1, 13.0);
  }

  #[test]
  fn test_combine_confidence_quadrature() {
    // c1: (8, 12) -> deviations of 2 from center 10
    // c2: (7, 13) -> deviations of 3 from center 10
    // Combined: sqrt(2² + 3²) = sqrt(13) ≈ 3.606
    let result = combine_confidence(10.0, (0.0, 20.0), Some((8.0, 12.0)), Some((7.0, 13.0)));
    let expected_dev = (4.0_f64 + 9.0).sqrt();
    assert_relative_eq!(result.0, 10.0 - expected_dev, epsilon = 1e-10);
    assert_relative_eq!(result.1, 10.0 + expected_dev, epsilon = 1e-10);
  }

  #[test]
  fn test_combine_confidence_clipped_to_limits() {
    // Large contributions that exceed limits
    let result = combine_confidence(10.0, (5.0, 15.0), Some((0.0, 20.0)), Some((0.0, 20.0)));
    // Quadrature would give sqrt(100 + 100) = 14.14 deviation
    // But limits clip to (5, 15)
    assert_relative_eq!(result.0, 5.0);
    assert_relative_eq!(result.1, 15.0);
  }
}
