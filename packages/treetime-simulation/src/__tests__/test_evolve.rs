#[cfg(test)]
mod tests {
  use crate::__tests__::test_support::tests::helpers::three_tip_tree;
  use crate::evolve::{SiteEvolution, evolve_sequences, evolve_sequences_with_rates, evolve_site, gamma_site_rates};
  use crate::substitution_model::SubstitutionModel;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_primitives::AlignmentRecord;
  use treetime_utils::assert_error;
  use treetime_utils::cast::ToF64Exact;
  use treetime_utils::pretty_assert_abs_diff_eq;
  use treetime_utils::sync::random::get_random_number_generator;

  #[rustfmt::skip]
  #[rstest]
  #[case::a(0)]
  #[case::c(1)]
  #[case::g(2)]
  #[case::t(3)]
  #[trace]
  fn test_evolve_site_without_rate_keeps_state(#[case] state: usize) {
    let model = SubstitutionModel::jc69(0.0).unwrap();
    let expected = SiteEvolution { state, jumps: 0 };
    let actual = evolve_site(state, 10.0, &model, &mut get_random_number_generator(4)).unwrap();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_evolve_sequences_without_rate_copies_root_sequence() {
    let tree = three_tip_tree();
    let model = SubstitutionModel::random(0.0, &mut get_random_number_generator(5)).unwrap();
    let durations = btreemap! { tree.root_inner => 1.0, tree.root_c => 2.0, tree.inner_a => 3.0, tree.inner_b => 4.0 };
    let simulated = evolve_sequences(
      &tree.graph,
      &durations,
      &model,
      200,
      &mut get_random_number_generator(6),
    )
    .unwrap();

    let root_sequence = &simulated.sequences[&tree.root];
    let expected = (
      vec![root_sequence; 5],
      btreemap! { tree.root_inner => 0, tree.root_c => 0, tree.inner_a => 0, tree.inner_b => 0 },
    );
    let actual = (
      simulated.sequences.values().collect::<Vec<_>>(),
      simulated.substitutions.clone(),
    );
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_evolve_sequences_with_unit_rates_reproduces_evolve_sequences() {
    let tree = three_tip_tree();
    let model = SubstitutionModel::jc69(1.0).unwrap();
    let durations = btreemap! { tree.root_inner => 0.1, tree.root_c => 0.2, tree.inner_a => 0.3, tree.inner_b => 0.4 };
    let plain = evolve_sequences(&tree.graph, &durations, &model, 50, &mut get_random_number_generator(8)).unwrap();
    let rated = evolve_sequences_with_rates(
      &tree.graph,
      &durations,
      &model,
      &[1.0; 50],
      &mut get_random_number_generator(8),
    )
    .unwrap();
    assert_eq!(plain, rated);
  }

  #[test]
  fn test_evolve_sequences_with_rates_keeps_sites_of_zero_rate() {
    let tree = three_tip_tree();
    let model = SubstitutionModel::jc69(1.0).unwrap();
    let durations = btreemap! { tree.root_inner => 5.0, tree.root_c => 5.0, tree.inner_a => 5.0, tree.inner_b => 5.0 };
    let rates = (0..100)
      .map(|site| if site % 2 == 0 { 0.0 } else { 1.0 })
      .collect::<Vec<_>>();
    let simulated = evolve_sequences_with_rates(
      &tree.graph,
      &durations,
      &model,
      &rates,
      &mut get_random_number_generator(9),
    )
    .unwrap();
    let root = &simulated.sequences[&tree.root];
    for sequence in simulated.sequences.values() {
      for site in (0..100).step_by(2) {
        assert_eq!(sequence[site], root[site], "site {site}");
      }
    }
  }

  #[test]
  fn test_gamma_site_rates_have_mean_one_and_shape_variance() {
    let alpha = 0.5;
    let rates = gamma_site_rates(alpha, 200_000, &mut get_random_number_generator(10)).unwrap();
    let n = rates.len().to_f64_exact();
    let mean = rates.iter().sum::<f64>() / n;
    let variance = rates.iter().map(|rate| (rate - mean).powi(2)).sum::<f64>() / n;
    pretty_assert_abs_diff_eq!(mean, 1.0, epsilon = 0.01);
    pretty_assert_abs_diff_eq!(variance, 1.0 / alpha, epsilon = 0.05);
  }

  #[test]
  fn test_evolve_sequences_leaf_records_name_leaf_sequences() {
    let tree = three_tip_tree();
    let model = SubstitutionModel::jc69(1.0).unwrap();
    let durations = btreemap! { tree.root_inner => 0.1, tree.root_c => 0.2, tree.inner_a => 0.3, tree.inner_b => 0.4 };
    let simulated = evolve_sequences(&tree.graph, &durations, &model, 20, &mut get_random_number_generator(7)).unwrap();

    let expected = vec![
      AlignmentRecord {
        name: "A".to_owned(),
        seq: simulated.sequences[&tree.tip_a].clone(),
      },
      AlignmentRecord {
        name: "B".to_owned(),
        seq: simulated.sequences[&tree.tip_b].clone(),
      },
      AlignmentRecord {
        name: "C".to_owned(),
        seq: simulated.sequences[&tree.tip_c].clone(),
      },
    ];
    let actual = simulated.leaf_records(&tree.graph, &tree.names);
    assert_eq!(expected, actual);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::invalid_state(     4, 1.0,           "State 4 is not one of the 4 nucleotide states")]
  #[case::negative_duration( 0, -1.0,          "Duration must be finite and non-negative, got -1")]
  #[case::nan_duration(      0, f64::NAN,      "Duration must be finite and non-negative, got NaN")]
  #[case::infinite_duration( 0, f64::INFINITY, "Duration must be finite and non-negative, got inf")]
  #[trace]
  fn test_evolve_site_rejects_invalid_input(#[case] state: usize, #[case] duration: f64, #[case] expected: &str) {
    let model = SubstitutionModel::jc69(1.0).unwrap();
    let result = evolve_site(state, duration, &model, &mut get_random_number_generator(0));
    assert_error!(result, expected);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::zero_length(      0, Some(1.0),           "Sequence length must be positive, got 0")]
  #[case::missing_duration( 10, None,               "Edge 3 has no duration")]
  #[case::negative_duration(10, Some(-1.0),         "Duration of edge 3 must be finite and non-negative, got -1")]
  #[case::nan_duration(     10, Some(f64::NAN),     "Duration of edge 3 must be finite and non-negative, got NaN")]
  #[case::infinite_duration(10, Some(f64::INFINITY), "Duration of edge 3 must be finite and non-negative, got inf")]
  #[trace]
  fn test_evolve_sequences_rejects_invalid_input(
    #[case] length: usize,
    #[case] last_duration: Option<f64>,
    #[case] expected: &str,
  ) {
    let tree = three_tip_tree();
    let model = SubstitutionModel::jc69(1.0).unwrap();
    let mut durations: BTreeMap<GraphEdgeKey, f64> =
      btreemap! { tree.root_inner => 0.1, tree.root_c => 0.2, tree.inner_a => 0.3 };
    if let Some(duration) = last_duration {
      durations.insert(tree.inner_b, duration);
    }
    let result = evolve_sequences(&tree.graph, &durations, &model, length, &mut get_random_number_generator(0));
    assert_error!(result, expected);
  }
}
