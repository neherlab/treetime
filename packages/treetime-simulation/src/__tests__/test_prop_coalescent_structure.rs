#[cfg(test)]
mod tests {
  use crate::coalescent::simulate_coalescent;
  use generators::{gen_sampling, gen_schedule};
  use itertools::Itertools;
  use proptest::{prop_assert_eq, proptest};
  use std::collections::BTreeMap;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_utils::sync::random::get_random_number_generator;

  proptest! {
    #[test]
    fn test_prop_coalescent_tips_match_sampling(
      schedule in gen_schedule(),
      sampling in gen_sampling(),
      seed in 0..u64::MAX,
    ) {
      let tree = simulate_coalescent(&schedule, &sampling, &mut get_random_number_generator(seed)).unwrap();

      let expected = sampling
        .iter()
        .flat_map(|batch| vec![batch.tau; batch.count])
        .sorted_by(f64::total_cmp)
        .collect_vec();
      let actual = tree
        .graph
        .get_leaves()
        .map(|leaf| tree.tau[&leaf.key()])
        .sorted_by(f64::total_cmp)
        .collect_vec();
      prop_assert_eq!(expected, actual);
    }

    #[test]
    fn test_prop_coalescent_tree_is_binary_with_one_root(
      schedule in gen_schedule(),
      sampling in gen_sampling(),
      seed in 0..u64::MAX,
    ) {
      let tree = simulate_coalescent(&schedule, &sampling, &mut get_random_number_generator(seed)).unwrap();

      let tip_count: usize = sampling.iter().map(|batch| batch.count).sum();
      let expected = (1, vec![2; tip_count - 1], vec![1; 2 * tip_count - 2]);
      let actual = (
        tree.graph.get_roots().count(),
        tree.graph.get_internal_nodes().map(|node| node.outbound().len()).collect_vec(),
        tree.graph.get_nodes().filter(|node| !node.is_root()).map(|node| node.inbound().len()).collect_vec(),
      );
      prop_assert_eq!(expected, actual);
    }

    #[test]
    fn test_prop_coalescent_parents_are_strictly_older_than_children(
      schedule in gen_schedule(),
      sampling in gen_sampling(),
      seed in 0..u64::MAX,
    ) {
      let tree = simulate_coalescent(&schedule, &sampling, &mut get_random_number_generator(seed)).unwrap();

      let durations = tree.edge_durations();
      let tau_differences: BTreeMap<GraphEdgeKey, f64> = tree
        .graph
        .get_edges()
        .map(|edge| (edge.key(), tree.tau[&edge.source()] - tree.tau[&edge.target()]))
        .collect();
      let invalid_durations = durations
        .iter()
        .filter(|(_, duration)| !(duration.is_finite() && **duration > 0.0))
        .map(|(edge, _)| *edge)
        .collect_vec();
      prop_assert_eq!((tau_differences, vec![]), (durations, invalid_durations));
    }
  }

  mod generators {
    use crate::coalescent::SampleBatch;
    use crate::tc_schedule::TcSchedule;
    use proptest::collection::vec;
    use proptest::prop_oneof;
    use proptest::strategy::{Just, Strategy};

    pub(super) fn gen_schedule() -> impl Strategy<Value = TcSchedule> {
      prop_oneof![
        Just(TcSchedule::constant(2.0).unwrap()),
        Just(TcSchedule::exponential(2.0, 1.5).unwrap()),
        Just(TcSchedule::piecewise(vec![1.0, 2.0], vec![2.0, 0.5, 1.5]).unwrap()),
      ]
    }

    pub(super) fn gen_sampling() -> impl Strategy<Value = Vec<SampleBatch>> {
      vec((0.0..3.0_f64, 1..5_usize), 1..5)
        .prop_map(|batches| {
          batches
            .into_iter()
            .map(|(tau, count)| SampleBatch { tau, count })
            .collect::<Vec<_>>()
        })
        .prop_filter("at least 2 tips", |batches| {
          batches.iter().map(|batch| batch.count).sum::<usize>() >= 2
        })
    }
  }
}
