#[cfg(test)]
mod tests {
  use crate::timetree::confidence::extract_confidence_intervals;
  use approx::assert_relative_eq;
  use helpers::add_named;
  use maplit::btreemap;
  use ndarray::Array1;
  use statrs::distribution::{ContinuousCDF, Normal};
  use std::collections::BTreeMap;
  use std::sync::Arc;
  use treetime_distribution::Distribution;
  use treetime_graph::graph::Graph;
  use treetime_utils::cast::ToF64Exact;

  #[test]
  fn test_extract_confidence_intervals_includes_unnamed_nodes() {
    let mut graph = Graph::new();
    let mut names = BTreeMap::new();
    let unnamed_key = add_named(&mut graph, &mut names, None);
    let named_key = add_named(&mut graph, &mut names, Some("named"));
    graph.build().unwrap();

    let state = helpers::state(
      &graph,
      &[(unnamed_key, Some(2020.0), None), (named_key, Some(2021.0), None)],
    );
    let intervals = extract_confidence_intervals(&graph, &state, &BTreeMap::new(), &names).unwrap();
    assert_eq!(2, intervals.len());
    assert_eq!("", intervals[0].name);
    assert_eq!("named", intervals[1].name);
  }

  #[test]
  fn test_extract_confidence_intervals_skips_nodes_without_time() {
    let mut graph = Graph::new();
    let mut names = BTreeMap::new();
    let no_time_key = add_named(&mut graph, &mut names, Some("no_time"));
    let has_time_key = add_named(&mut graph, &mut names, Some("has_time"));
    graph.build().unwrap();

    let state = helpers::state(&graph, &[(no_time_key, None, None), (has_time_key, Some(2021.0), None)]);
    let intervals = extract_confidence_intervals(&graph, &state, &BTreeMap::new(), &names).unwrap();
    assert_eq!(1, intervals.len());
    assert_eq!("has_time", intervals[0].name);
  }

  #[test]
  fn test_extract_confidence_intervals_uses_date_as_fallback() {
    let mut graph = Graph::new();
    let mut names = BTreeMap::new();
    let key = add_named(&mut graph, &mut names, Some("node_a"));
    graph.build().unwrap();

    let state = helpers::state(&graph, &[(key, Some(2020.5), None)]);
    let intervals = extract_confidence_intervals(&graph, &state, &BTreeMap::new(), &names).unwrap();
    assert_eq!(1, intervals.len());
    assert_relative_eq!(intervals[0].date, 2020.5);
    assert_relative_eq!(intervals[0].lower, 2020.5);
    assert_relative_eq!(intervals[0].upper, 2020.5);
  }

  #[test]
  fn test_extract_confidence_intervals_flat_distribution_extends_from_latest_date() {
    let mut graph = Graph::new();
    let mut names = BTreeMap::new();
    let dist = Arc::new(Distribution::range((2019.0, 2021.0), 0.0));
    let key = add_named(&mut graph, &mut names, Some("node_a"));
    graph.build().unwrap();

    let state = helpers::state(&graph, &[(key, Some(2020.0), Some(dist))]);
    let intervals = extract_confidence_intervals(&graph, &state, &BTreeMap::new(), &names).unwrap();
    assert_eq!(1, intervals.len());
    assert_relative_eq!(intervals[0].date, 2020.0);
    let uniform_mass_from_latest_date = (2021.0 - 0.9 * 2.0, 2021.0);
    assert_relative_eq!(intervals[0].lower, uniform_mass_from_latest_date.0, epsilon = 1e-12);
    assert_relative_eq!(intervals[0].upper, uniform_mass_from_latest_date.1, epsilon = 1e-12);
  }

  #[test]
  fn test_extract_confidence_intervals_sorted_by_key() {
    let mut graph = Graph::new();
    let mut names = BTreeMap::new();
    let zebra_key = add_named(&mut graph, &mut names, Some("zebra"));
    let alpha_key = add_named(&mut graph, &mut names, Some("alpha"));
    let middle_key = add_named(&mut graph, &mut names, Some("middle"));
    graph.build().unwrap();

    let state = helpers::state(
      &graph,
      &[
        (zebra_key, Some(2020.0), None),
        (alpha_key, Some(2021.0), None),
        (middle_key, Some(2022.0), None),
      ],
    );
    let intervals = extract_confidence_intervals(&graph, &state, &BTreeMap::new(), &names).unwrap();
    assert_eq!(3, intervals.len());
    assert_eq!("zebra", intervals[0].name);
    assert_eq!("alpha", intervals[1].name);
    assert_eq!("middle", intervals[2].name);
  }

  #[test]
  fn test_extract_confidence_intervals_rate_only() {
    let mut graph = Graph::new();
    let mut names = BTreeMap::new();
    let key = add_named(&mut graph, &mut names, Some("node_a"));
    graph.build().unwrap();
    let rate_map = btreemap! { key => [2009.0, 2010.0, 2011.0] };

    let state = helpers::state(&graph, &[(key, Some(2010.0), None)]);
    let intervals = extract_confidence_intervals(&graph, &state, &rate_map, &names).unwrap();
    assert_eq!(1, intervals.len());
    assert_relative_eq!(intervals[0].date, 2010.0);
    assert_relative_eq!(intervals[0].lower, 2008.355146, epsilon = 1e-4);
    assert_relative_eq!(intervals[0].upper, 2011.644854, epsilon = 1e-4);
  }

  #[test]
  fn test_extract_confidence_intervals_combined_wider_than_either() {
    let mut graph = Graph::new();
    let mut names = BTreeMap::new();
    let (mean, std_dev): (f64, f64) = (2010.0, 0.5);
    let t = Array1::linspace(mean - 8.0 * std_dev, mean + 8.0 * std_dev, 16001);
    let y = t.mapv(|x| (x - mean).powi(2) / (2.0 * std_dev * std_dev));
    let dist = Arc::new(Distribution::function(t, y).unwrap());
    let key = add_named(&mut graph, &mut names, Some("node_a"));
    graph.build().unwrap();
    let rate_map = btreemap! { key => [2009.0, 2010.0, 2011.0] };

    let state = helpers::state(&graph, &[(key, Some(mean), Some(dist))]);
    let intervals = extract_confidence_intervals(&graph, &state, &rate_map, &names).unwrap();
    assert_eq!(1, intervals.len());

    let z = Normal::standard().inverse_cdf(0.95);
    let normal_hpd_half_width = z * std_dev;
    let rate_half_width = z * 1.0;
    let half_width = normal_hpd_half_width.hypot(rate_half_width);
    assert_relative_eq!(intervals[0].lower, mean - half_width, epsilon = 1e-7);
    assert_relative_eq!(intervals[0].upper, mean + half_width, epsilon = 1e-7);
  }

  #[test]
  fn test_extract_confidence_intervals_clamps_when_date_outside_rate_ci() {
    let mut graph = Graph::new();
    let mut names = BTreeMap::new();
    let key = add_named(&mut graph, &mut names, Some("node_a"));
    graph.build().unwrap();
    let rate_map = btreemap! { key => [2019.9, 2020.0, 2020.1] };

    let state = helpers::state(&graph, &[(key, Some(2020.5), None)]);
    let intervals = extract_confidence_intervals(&graph, &state, &rate_map, &names).unwrap();
    assert_eq!(1, intervals.len());
    assert!(intervals[0].lower <= intervals[0].date);
    assert!(intervals[0].date <= intervals[0].upper);
    assert_relative_eq!(intervals[0].upper, 2020.5);
    assert_relative_eq!(intervals[0].lower, 2019.8355, epsilon = 1e-3);
  }

  #[test]
  fn test_extract_confidence_intervals_clamps_when_date_below_rate_ci() {
    let mut graph = Graph::new();
    let mut names = BTreeMap::new();
    let key = add_named(&mut graph, &mut names, Some("node_a"));
    graph.build().unwrap();
    let rate_map = btreemap! { key => [2019.9, 2020.0, 2020.1] };

    let state = helpers::state(&graph, &[(key, Some(2019.5), None)]);
    let intervals = extract_confidence_intervals(&graph, &state, &rate_map, &names).unwrap();
    assert_eq!(1, intervals.len());
    assert!(intervals[0].lower <= intervals[0].date);
    assert!(intervals[0].date <= intervals[0].upper);
    assert_relative_eq!(intervals[0].lower, 2019.5);
    assert_relative_eq!(intervals[0].upper, 2020.1645, epsilon = 1e-3);
  }

  #[test]
  fn test_extract_confidence_intervals_skewed_distribution_hpd() {
    let n_points = 10001;
    let x_min = 0.0;
    let dx = 10.0 / (n_points.to_f64_exact() - 1.0);
    let y = Array1::from_shape_fn(n_points, |i| x_min + i.to_f64_exact() * dx);

    let dist_fn = treetime_distribution::DistributionFunction::from_start_dx_values(x_min, dx, y).unwrap();
    let dist = Distribution::Function(dist_fn);
    let peak_time = dist.likely_time().unwrap().unwrap();

    let mut graph = Graph::new();
    let mut names = BTreeMap::new();
    let key = add_named(&mut graph, &mut names, Some("skewed"));
    graph.build().unwrap();

    let state = helpers::state(&graph, &[(key, Some(peak_time), Some(Arc::new(dist)))]);
    let intervals = extract_confidence_intervals(&graph, &state, &BTreeMap::new(), &names).unwrap();
    assert_eq!(1, intervals.len());

    let hpd_lower = 0.0;
    let truncated_exponential_quantile = -(1.0 - 0.9 * (1.0 - (-10.0_f64).exp())).ln();
    assert_relative_eq!(intervals[0].lower, hpd_lower, epsilon = 1e-12);
    assert_relative_eq!(intervals[0].upper, truncated_exponential_quantile, epsilon = 1e-7);
  }

  mod helpers {
    use crate::timetree::inference::result::NodePosterior;
    use std::collections::BTreeMap;
    use std::sync::Arc;
    use treetime_distribution::{Distribution, NegLog};
    use treetime_graph::graph::Graph;
    use treetime_graph::node::GraphNodeKey;

    pub(super) type NodeTimeEntry = (GraphNodeKey, Option<f64>, Option<Arc<Distribution<NegLog>>>);

    pub(super) fn state(graph: &Graph, entries: &[NodeTimeEntry]) -> BTreeMap<GraphNodeKey, NodePosterior> {
      let mut posterior: BTreeMap<GraphNodeKey, NodePosterior> = graph
        .get_nodes()
        .map(|node| (node.key(), NodePosterior::default()))
        .collect();
      for (key, time, distribution) in entries {
        let likely_time = distribution.as_deref().and_then(|distribution| {
          distribution
            .likely_time()
            .expect("fixture distribution has a likely time")
        });
        let node = NodePosterior {
          distribution: distribution.clone(),
          likely_time,
          time: *time,
          contradicted: false,
        };
        posterior.insert(*key, node);
      }
      posterior
    }

    pub(super) fn add_named(
      graph: &mut Graph,
      names: &mut BTreeMap<GraphNodeKey, Option<String>>,
      name: Option<&str>,
    ) -> GraphNodeKey {
      let key = graph.add_node();
      names.insert(key, name.map(|n| n.to_owned()));
      key
    }
  }
}
