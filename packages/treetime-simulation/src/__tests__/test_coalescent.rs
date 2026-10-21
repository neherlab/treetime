#[cfg(test)]
mod tests {
  use crate::coalescent::{PRESENT, SampleBatch, simulate_coalescent};
  use crate::tc_schedule::TcSchedule;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::assert_error;
  use treetime_utils::sync::random::get_random_number_generator;

  #[test]
  fn test_coalescent_names_tips_internal_nodes_and_root() {
    let schedule = TcSchedule::constant(1.0).unwrap();
    let sampling = [SampleBatch { tau: 0.0, count: 2 }, SampleBatch { tau: 0.5, count: 1 }];
    let tree = simulate_coalescent(&schedule, &sampling, &mut get_random_number_generator(1)).unwrap();

    let root = tree.graph.root_key().unwrap();
    let expected = (
      "node_root",
      vec!["node_0001", "node_root", "tip_0001", "tip_0002", "tip_0003"],
    );
    let actual = (
      &tree.names[root],
      tree
        .graph
        .get_nodes()
        .map(|node| &tree.names[node.key()])
        .sorted()
        .collect_vec(),
    );
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_coalescent_dates_count_back_from_present() {
    let schedule = TcSchedule::constant(1.0).unwrap();
    let sampling = [SampleBatch { tau: 0.0, count: 1 }, SampleBatch { tau: 1.5, count: 1 }];
    let tree = simulate_coalescent(&schedule, &sampling, &mut get_random_number_generator(2)).unwrap();

    let expected = vec![PRESENT, PRESENT - 1.5];
    let actual = tree.graph.get_leaves().map(|leaf| tree.date(leaf.key())).collect_vec();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_coalescent_sorts_batches_by_time() {
    let schedule = TcSchedule::constant(1.0).unwrap();
    let sorted = [SampleBatch { tau: 0.0, count: 3 }, SampleBatch { tau: 1.0, count: 2 }];
    let unsorted = [SampleBatch { tau: 1.0, count: 2 }, SampleBatch { tau: 0.0, count: 3 }];

    let expected = simulate_coalescent(&schedule, &sorted, &mut get_random_number_generator(3)).unwrap();
    let actual = simulate_coalescent(&schedule, &unsorted, &mut get_random_number_generator(3)).unwrap();
    assert_eq!(expected.tau, actual.tau);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::no_batch(          vec![],                                                                    "Sampling must have at least one batch")]
  #[case::batch_without_tips(vec![SampleBatch { tau: 0.0, count: 3 }, SampleBatch { tau: 1.0, count: 0 }], "Sample batch 1 at time 1 has no tips")]
  #[case::single_tip(        vec![SampleBatch { tau: 0.0, count: 1 }],                                  "Sampling must have at least 2 tips in total, got 1")]
  #[case::negative_time(     vec![SampleBatch { tau: -1.0, count: 2 }],                                 "Sample batch 0 has time -1, expected a finite non-negative number of years before the most recent sample")]
  #[case::nan_time(          vec![SampleBatch { tau: 0.0, count: 2 }, SampleBatch { tau: f64::NAN, count: 2 }], "Sample batch 1 has time NaN, expected a finite non-negative number of years before the most recent sample")]
  #[case::infinite_time(     vec![SampleBatch { tau: f64::INFINITY, count: 2 }],                        "Sample batch 0 has time inf, expected a finite non-negative number of years before the most recent sample")]
  #[trace]
  fn test_coalescent_rejects_invalid_sampling(#[case] sampling: Vec<SampleBatch>, #[case] expected: &str) {
    let schedule = TcSchedule::constant(1.0).unwrap();
    let result = simulate_coalescent(&schedule, &sampling, &mut get_random_number_generator(0));
    assert_error!(result, expected);
  }
}
